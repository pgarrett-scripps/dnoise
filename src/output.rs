//! Path checks and owned temporary output. Only completed runs replace files.

use crate::{DnoiseError, Result};
use std::fs::{self, OpenOptions};
use std::path::{Component, Path, PathBuf};
use tempfile::{Builder, TempDir};

/// Resolve existing ancestors (including symlinks) and normalize missing suffixes.
pub fn resolved_path(path: &Path) -> Result<PathBuf> {
    let absolute = if path.is_absolute() {
        path.to_owned()
    } else {
        std::env::current_dir()?.join(path)
    };
    let mut result = PathBuf::new();
    for part in absolute.components() {
        match part {
            // A Windows drive/UNC prefix is not a filesystem path on its own,
            // especially the verbatim prefix returned by canonicalize(). Wait
            // for RootDir before querying it.
            Component::Prefix(prefix) => result.push(prefix.as_os_str()),
            Component::ParentDir => {
                result.pop();
            }
            Component::CurDir => {}
            component => {
                result.push(component.as_os_str());
                match fs::symlink_metadata(&result) {
                    Ok(_) => result = result.canonicalize()?,
                    Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
                    Err(e) => return Err(e.into()),
                }
            }
        }
    }
    Ok(result)
}

/// Reject equal paths and either directory containing the other.
pub fn check_disjoint(input: &Path, output: &Path) -> Result<()> {
    let a = resolved_path(input)?;
    let b = resolved_path(output)?;
    if a.starts_with(&b) || b.starts_with(&a) {
        return Err(DnoiseError::InvalidInput(format!(
            "input and output paths overlap: {} and {}",
            input.display(),
            output.display()
        )));
    }
    Ok(())
}

pub(crate) struct OutputTransaction {
    staging: TempDir,
    destination: PathBuf,
    lock: PathBuf,
    overwrite: bool,
}

impl OutputTransaction {
    pub(crate) fn begin(
        input: &Path,
        output: &Path,
        overwrite: bool,
        in_place: bool,
    ) -> Result<Self> {
        if !in_place {
            check_disjoint(input, output)?;
        }
        let destination = resolved_path(output)?;
        if destination.exists() && !overwrite {
            return Err(DnoiseError::OutputExists(destination));
        }
        if destination.exists() && !destination.is_dir() {
            return Err(DnoiseError::InvalidInput(
                "output must be a directory".into(),
            ));
        }
        let parent = destination
            .parent()
            .ok_or_else(|| DnoiseError::InvalidInput("output has no parent".into()))?;
        // Do not create arbitrary ancestor directories during a failed run.
        if !parent.is_dir() {
            return Err(DnoiseError::InvalidInput(format!(
                "output parent does not exist: {}",
                parent.display()
            )));
        }
        let mut lock_name = destination
            .file_name()
            .ok_or_else(|| DnoiseError::InvalidInput("output has no name".into()))?
            .to_os_string();
        lock_name.push(".dnoise-lock");
        let lock = parent.join(lock_name);
        let lock_file = OpenOptions::new().write(true).create_new(true).open(&lock)
            .map_err(|e| DnoiseError::InvalidInput(format!("cannot acquire output lock {}: {e}; if a previous run crashed, inspect its temporary/backup folders before removing the lock", lock.display())))?;
        drop(lock_file);
        let staging = match Builder::new().prefix(".dnoise-stage-").tempdir_in(parent) {
            Ok(dir) => dir,
            Err(e) => {
                let _ = fs::remove_file(&lock);
                return Err(e.into());
            }
        };
        Ok(Self {
            staging,
            destination,
            lock,
            overwrite,
        })
    }

    pub(crate) fn path(&self) -> &Path {
        self.staging.path()
    }

    pub(crate) fn commit(self) -> Result<()> {
        self.commit_with(|a, b| fs::rename(a, b))
    }

    fn commit_with(
        self,
        mut rename: impl FnMut(&Path, &Path) -> std::io::Result<()>,
    ) -> Result<()> {
        let mut backup = None;
        if self.destination.exists() {
            if !self.overwrite {
                return Err(DnoiseError::OutputExists(self.destination.clone()));
            }
            let dir = Builder::new()
                .prefix(".dnoise-backup-")
                .tempdir_in(self.destination.parent().unwrap())?;
            rename(&self.destination, &dir.path().join("original"))?;
            // From this point on Drop must never delete the original on errors.
            backup = Some(dir.keep());
        }
        if let Err(error) = rename(self.path(), &self.destination) {
            if let Some(backup) = &backup {
                if let Err(restore) = rename(&backup.join("original"), &self.destination) {
                    return Err(DnoiseError::Recovery(format!(
                        "install failed: {error}; restore failed: {restore}; original preserved at {}",
                        backup.join("original").display()
                    )));
                }
                let _ = fs::remove_dir(backup);
            }
            let recovery = if backup.is_some() {
                "existing output restored"
            } else {
                "no existing output was replaced"
            };
            return Err(DnoiseError::Recovery(format!(
                "install failed: {error}; {recovery}"
            )));
        }
        if let Some(backup) = backup
            && let Err(e) = fs::remove_dir_all(&backup)
        {
            tracing::warn!(path = %backup.display(), "output installed; backup cleanup failed: {e}");
        }
        Ok(())
    }
}

impl Drop for OutputTransaction {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.lock);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn setup() -> (TempDir, PathBuf, PathBuf) {
        let root = tempfile::tempdir().unwrap();
        let input = root.path().join("input.d");
        let output = root.path().join("output.d");
        fs::create_dir(&input).unwrap();
        fs::create_dir(&output).unwrap();
        fs::write(output.join("original"), b"keep").unwrap();
        (root, input, output)
    }
    #[test]
    fn rejects_equal_nested_and_ancestor_paths() {
        let (_root, input, _) = setup();
        for output in [
            input.clone(),
            input.join("nested.d"),
            input.parent().unwrap().to_owned(),
        ] {
            assert!(check_disjoint(&input, &output).is_err());
        }
    }
    #[test]
    fn resolves_canonical_paths_and_missing_children() {
        let (_root, input, _) = setup();
        let canonical = input.canonicalize().unwrap();
        assert_eq!(resolved_path(&canonical).unwrap(), canonical);
        assert_eq!(
            resolved_path(&canonical.join("missing.d")).unwrap(),
            canonical.join("missing.d")
        );
        assert!(check_disjoint(&input, &canonical.join("missing.d")).is_err());
    }
    #[test]
    fn failed_run_preserves_output_and_releases_lock() {
        let (_root, input, output) = setup();
        let tx = OutputTransaction::begin(&input, &output, true, false).unwrap();
        fs::write(tx.path().join("partial"), b"partial").unwrap();
        assert!(OutputTransaction::begin(&input, &output, true, false).is_err());
        let staging = tx.path().to_owned();
        drop(tx);
        assert!(!staging.exists());
        assert_eq!(fs::read(output.join("original")).unwrap(), b"keep");
        assert!(OutputTransaction::begin(&input, &output, true, false).is_ok());
    }
    #[test]
    fn failed_install_restores_original() {
        let (_root, input, output) = setup();
        let tx = OutputTransaction::begin(&input, &output, true, false).unwrap();
        let stage = tx.path().to_owned();
        let result = tx.commit_with(|a, b| {
            if a == stage {
                Err(std::io::Error::other("injected install failure"))
            } else {
                fs::rename(a, b)
            }
        });
        assert!(result.is_err());
        assert_eq!(fs::read(output.join("original")).unwrap(), b"keep");
    }
    #[test]
    fn failed_restore_keeps_backup_and_reports_location() {
        let (root, input, output) = setup();
        let tx = OutputTransaction::begin(&input, &output, true, false).unwrap();
        let mut calls = 0;
        let error = tx
            .commit_with(|a, b| {
                calls += 1;
                if calls > 1 {
                    Err(std::io::Error::other("injected failure"))
                } else {
                    fs::rename(a, b)
                }
            })
            .unwrap_err();
        assert!(error.to_string().contains("original preserved at"));
        let backup = fs::read_dir(root.path())
            .unwrap()
            .flatten()
            .find(|e| {
                e.file_name()
                    .to_string_lossy()
                    .starts_with(".dnoise-backup-")
            })
            .unwrap();
        assert_eq!(
            fs::read(backup.path().join("original/original")).unwrap(),
            b"keep"
        );
    }
    #[cfg(unix)]
    #[test]
    fn rejects_symlink_alias_and_normalizes_missing_suffix() {
        let (root, input, _) = setup();
        let alias = root.path().join("alias");
        std::os::unix::fs::symlink(&input, &alias).unwrap();
        assert!(check_disjoint(&input, &alias).is_err());
        assert!(check_disjoint(&input, &alias.join("new.d")).is_err());
    }
}
