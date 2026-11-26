//! Path checks and owned temporary output. Only completed runs replace files.

use crate::{DnoiseError, Result};
use same_file::Handle;
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

/// Reject equal paths and either directory containing the other, including the
/// case where the two paths reach one directory by different routes.
///
/// Comparing paths is not enough on its own. `canonicalize` resolves symlinks
/// but not bind mounts, hard links or drive substitutions, so two paths can
/// normalize to different strings and still name the same directory. Replacing
/// such an output deletes the input, so the identity of the files themselves is
/// compared as well.
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
    check_not_aliased(&a, &b)
}

/// Identity of an existing file or directory, independent of the path used to
/// reach it: device and inode on Unix, volume serial number and file index on
/// Windows. `None` when the path does not exist.
fn identity(path: &Path) -> Result<Option<Handle>> {
    match Handle::from_path(path) {
        Ok(handle) => Ok(Some(handle)),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
        // A path we may not open cannot be compared. Treat it as unknown rather
        // than failing a run that never touches it; the open of the real input
        // or output reports the permission error with its own path.
        Err(e) if e.kind() == std::io::ErrorKind::PermissionDenied => Ok(None),
        Err(e) => Err(e.into()),
    }
}

/// Reject any pair where one path and the other, or one of the other's existing
/// ancestors, are the same file on disk. Both directions are checked: the
/// output can be the input reached another way, and it can also sit inside it.
fn check_not_aliased(input: &Path, output: &Path) -> Result<()> {
    for (one, two) in [(input, output), (output, input)] {
        let Some(id) = identity(one)? else { continue };
        for ancestor in two.ancestors() {
            if let Some(other) = identity(ancestor)?
                && other == id
            {
                return Err(DnoiseError::InvalidInput(format!(
                    "input and output reach the same file on disk by two different \
                     paths, so writing the output would destroy the input: {} and {} \
                     (a bind mount, a hard link or a substituted drive)",
                    one.display(),
                    ancestor.display()
                )));
            }
        }
    }
    Ok(())
}

/// A destination may only be replaced when it is empty or is itself a `.d`
/// folder. `--force` used to move any existing directory aside and then delete
/// it, so an output path that pointed at an ordinary folder destroyed it.
pub(crate) fn check_replaceable(destination: &Path) -> Result<()> {
    if fs::read_dir(destination)?.next().is_none() || destination.join("analysis.tdf").is_file() {
        return Ok(());
    }
    Err(DnoiseError::InvalidInput(format!(
        "refusing to replace {}: it is not empty and holds no analysis.tdf, so it is \
         not a .d folder left by an earlier run; point the output somewhere else, or \
         remove that folder yourself first",
        destination.display()
    )))
}

pub(crate) struct OutputTransaction {
    staging: TempDir,
    destination: PathBuf,
    lock: PathBuf,
    overwrite: bool,
    /// The resolved input, kept so the destructive step in `commit` can repeat
    /// the checks against whatever the destination resolves to by then. `None`
    /// for an in-place run, where input and output are the same folder on
    /// purpose.
    guarded_input: Option<PathBuf>,
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
        if destination.exists() {
            check_replaceable(&destination)?;
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
            guarded_input: if in_place {
                None
            } else {
                Some(resolved_path(input)?)
            },
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
            // The run took minutes, so repeat both checks at the one moment
            // that is destructive: the destination may have been replaced by an
            // alias of the input, or by a folder that is not a .d, since begin.
            if let Some(input) = &self.guarded_input {
                check_not_aliased(input, &self.destination)?;
            }
            check_replaceable(&self.destination)?;
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
        // A replaceable destination is empty or holds analysis.tdf. These
        // fixtures stand in for a previous run's output.
        fs::write(output.join("analysis.tdf"), b"tdf").unwrap();
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
    /// `--force` once moved any existing directory aside and deleted it, so an
    /// output path that pointed at an ordinary folder destroyed it.
    #[test]
    fn force_refuses_a_destination_that_is_not_a_d_folder() {
        let (root, input, _) = setup();
        let destination = root.path().join("notes");
        fs::create_dir(&destination).unwrap();
        fs::write(destination.join("thesis.txt"), b"irreplaceable").unwrap();
        let error = OutputTransaction::begin(&input, &destination, true, false)
            .err()
            .unwrap();
        assert!(
            error.to_string().contains("holds no analysis.tdf"),
            "{error}"
        );
        assert_eq!(
            fs::read(destination.join("thesis.txt")).unwrap(),
            b"irreplaceable"
        );
        // An empty folder and a previous .d output are both replaceable.
        fs::remove_file(destination.join("thesis.txt")).unwrap();
        assert!(OutputTransaction::begin(&input, &destination, true, false).is_ok());
        fs::write(destination.join("analysis.tdf"), b"tdf").unwrap();
        fs::write(destination.join("analysis.tdf_bin"), b"bin").unwrap();
        assert!(OutputTransaction::begin(&input, &destination, true, false).is_ok());
    }

    /// The destination can become unreplaceable while the run is in progress,
    /// and commit is the step that deletes it.
    #[test]
    fn commit_refuses_a_destination_that_stopped_being_a_d_folder() {
        let (_root, input, output) = setup();
        let tx = OutputTransaction::begin(&input, &output, true, false).unwrap();
        fs::remove_file(output.join("analysis.tdf")).unwrap();
        let error = tx.commit().unwrap_err();
        assert!(
            error.to_string().contains("holds no analysis.tdf"),
            "{error}"
        );
        assert_eq!(fs::read(output.join("original")).unwrap(), b"keep");
    }

    /// Two paths can reach one directory through a bind mount or a substituted
    /// drive and still canonicalize differently, so the path comparison passes
    /// and the install deletes the input. A directory bind mount needs root, so
    /// the always-on case here is a hard-linked file, which has the same shape:
    /// one inode, two paths, neither a prefix of the other.
    #[cfg(unix)]
    #[test]
    fn rejects_same_file_reached_by_two_paths() {
        let (root, _, _) = setup();
        let one = root.path().join("one.json");
        let two = root.path().join("two.json");
        fs::write(&one, b"{}").unwrap();
        fs::hard_link(&one, &two).unwrap();
        assert_ne!(one.canonicalize().unwrap(), two.canonicalize().unwrap());
        for (a, b) in [(&one, &two), (&two, &one)] {
            let error = check_disjoint(a, b).unwrap_err();
            assert!(
                error.to_string().contains("the same file on disk"),
                "{error}"
            );
        }
        // A distinct file with the same contents is not an alias.
        let other = root.path().join("other.json");
        fs::write(&other, b"{}").unwrap();
        assert!(check_disjoint(&one, &other).is_ok());
    }

    #[cfg(unix)]
    #[test]
    fn identity_is_none_for_a_missing_path_and_equal_through_an_alias() {
        let (root, input, _) = setup();
        let one = root.path().join("one.json");
        let two = root.path().join("two.json");
        fs::write(&one, b"{}").unwrap();
        fs::hard_link(&one, &two).unwrap();
        assert_eq!(identity(&one).unwrap(), identity(&two).unwrap());
        assert_ne!(identity(&one).unwrap(), identity(&input).unwrap());
        assert!(identity(&root.path().join("missing")).unwrap().is_none());
    }

    /// Run against a real directory alias (the one on this machine is
    /// `/mnt/data/Data` and `~/Data`):
    /// `DNOISE_ALIAS_A=/mnt/data/Data DNOISE_ALIAS_B=$HOME/Data cargo test -- --ignored`
    #[test]
    #[ignore = "needs two paths that are the same directory, e.g. a bind mount"]
    fn rejects_a_real_directory_alias() {
        let (a, b) = match (
            std::env::var_os("DNOISE_ALIAS_A"),
            std::env::var_os("DNOISE_ALIAS_B"),
        ) {
            (Some(a), Some(b)) => (PathBuf::from(a), PathBuf::from(b)),
            _ => panic!("set DNOISE_ALIAS_A and DNOISE_ALIAS_B to two paths of one directory"),
        };
        assert_ne!(
            a.canonicalize().unwrap(),
            b.canonicalize().unwrap(),
            "the two paths canonicalize the same, so this is not the aliasing case"
        );
        // The same directory as input and as output.
        let error = check_disjoint(&a, &b).unwrap_err();
        assert!(
            error.to_string().contains("the same file on disk"),
            "{error}"
        );
        // An output that does not exist yet, under an aliased ancestor.
        let nested = check_disjoint(&a, &b.join("output.d")).unwrap_err();
        assert!(
            nested.to_string().contains("the same file on disk"),
            "{nested}"
        );
        // An input under one path, the output folder holding it under the other.
        let parent = check_disjoint(&a.join("run.d"), &b).unwrap_err();
        assert!(
            parent.to_string().contains("the same file on disk"),
            "{parent}"
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
