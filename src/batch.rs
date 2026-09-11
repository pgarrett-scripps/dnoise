//! Portable batch manifests shared by the CLI and desktop queue.
use crate::{DnoiseError, Result, config::Config};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

/// One independently configured input/output pair.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Job {
    /// Input acquisition directory.
    pub input: PathBuf,
    /// Output directory (must not overlap any batch input/output).
    pub output: PathBuf,
    /// Processing recipe; omitted fields use normal CLI defaults.
    #[serde(default)]
    pub config: Config,
}

/// Versioned JSON batch manifest. Relative paths are relative to the manifest.
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Manifest {
    /// Schema version (currently 1).
    pub schema_version: u32,
    /// Ordered jobs; execution continues past per-file processing errors.
    pub jobs: Vec<Job>,
}

impl Manifest {
    /// Load a manifest and resolve paths independently of the process directory.
    pub fn load(path: &Path) -> Result<Self> {
        let mut manifest: Self = serde_json::from_reader(std::fs::File::open(path)?)
            .map_err(|e| DnoiseError::InvalidInput(e.to_string()))?;
        if manifest.schema_version != 1 || manifest.jobs.is_empty() {
            return Err(DnoiseError::InvalidInput(
                "expected schema_version 1 and at least one job".into(),
            ));
        }
        let parent = path.canonicalize()?.parent().unwrap().to_owned();
        for job in &mut manifest.jobs {
            if job.input.is_relative() {
                job.input = parent.join(&job.input);
            }
            if job.output.is_relative() {
                job.output = parent.join(&job.output);
            }
        }
        Ok(manifest)
    }
}

/// Preflight the entire batch before any job writes, including cross-job collisions.
pub fn preflight(jobs: &[Job], force: bool) -> Result<()> {
    for (i, job) in jobs.iter().enumerate() {
        crate::validation::inspect(&job.input, false)?;
        job.config.resolve(&job.input)?;
        let destination = crate::output::resolved_path(&job.output)?;
        if destination.exists() && !force {
            return Err(DnoiseError::OutputExists(destination));
        }
        if destination.exists() && !destination.is_dir() {
            return Err(DnoiseError::InvalidInput(
                "output must be a directory".into(),
            ));
        }
        if !destination.parent().is_some_and(Path::is_dir) {
            return Err(DnoiseError::InvalidInput(format!(
                "output parent does not exist: {}",
                job.output.display()
            )));
        }
        for other in jobs {
            crate::output::check_disjoint(&other.input, &job.output)?;
        }
        for other in &jobs[..i] {
            crate::output::check_disjoint(&other.output, &job.output)?;
        }
    }
    Ok(())
}
