//! Inspection commands; ordinary denoising retains its positional CLI.
use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "dnoise")]
struct Commands {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Execute a JSON batch manifest, continuing past per-file errors.
    Batch {
        manifest: PathBuf,
        /// Replace existing outputs only after successful processing.
        #[arg(long)]
        force: bool,
        /// Estimate without writing acquisitions.
        #[arg(long)]
        dry_run: bool,
        /// Skip full input/output decoding checks for every job; keep structural checks.
        #[arg(long, alias = "no-validation", conflicts_with = "validate")]
        skip_validation: bool,
        /// Enable full validation for every job, overriding manifest settings.
        #[arg(long, conflicts_with = "skip_validation")]
        validate: bool,
        /// Save results JSON (must not exist); results are also printed on stdout.
        #[arg(long)]
        report: Option<PathBuf>,
    },
    /// Validate SQLite metadata and every decoded type-2 frame without writing.
    Validate { input: PathBuf },
    /// Print dnoise processing history; absence does not prove the input is raw.
    Metadata { input: PathBuf },
}

pub fn dispatch() -> Option<anyhow::Result<()>> {
    let name = std::env::args_os().nth(1)?;
    if name != "validate" && name != "metadata" && name != "batch" {
        return None;
    }
    Some((|| {
        match Commands::parse().command {
            Command::Batch {
                manifest,
                force,
                dry_run,
                skip_validation,
                validate,
                report,
            } => {
                let mut batch = dnoise::batch::Manifest::load(&manifest)?;
                for job in &mut batch.jobs {
                    if skip_validation {
                        job.config.skip_validation = Some(true);
                    } else if validate {
                        job.config.skip_validation = Some(false);
                    }
                }
                if !dry_run {
                    dnoise::batch::preflight(&batch.jobs, force)?;
                }
                if let Some(path) = &report {
                    dnoise::output::check_disjoint(&manifest, path)?;
                    for job in &batch.jobs {
                        dnoise::output::check_disjoint(&job.input, path)?;
                        dnoise::output::check_disjoint(&job.output, path)?;
                    }
                    if path.exists() {
                        anyhow::bail!("report already exists: {}", path.display());
                    }
                }
                let mut results = Vec::new();
                let mut failures = 0;
                for job in &batch.jobs {
                    let result = (|| -> anyhow::Result<serde_json::Value> {
                        let resolved = job.config.resolve(&job.input)?;
                        let stages = resolved.stages();
                        let options = dnoise::RunOptions {
                            force,
                            dry_run,
                            crop: (!resolved.crop.is_empty()).then_some(&resolved.crop),
                            crop_only: resolved.crop_only,
                            frame_batch_size: job.config.frame_batch_size,
                            skip_validation: job.config.skip_validation.unwrap_or(false),
                            ..Default::default()
                        };
                        let stats = job.config.in_thread_pool(|| {
                            dnoise::denoise_with_options(
                                &job.input,
                                &job.output,
                                &resolved.filter,
                                &stages,
                                &options,
                                |_| {},
                            )
                        })?;
                        Ok(dnoise::provenance::report(
                            &job.input,
                            &resolved.filter,
                            &stages,
                            &options,
                            &stats,
                        ))
                    })();
                    match result {
                        Ok(run)=>results.push(serde_json::json!({"input":job.input,"output":job.output,"success":true,"run":run})),
                        Err(e)=>{failures+=1;results.push(serde_json::json!({"input":job.input,"output":job.output,"success":false,"error":e.to_string()}));}
                    }
                }
                let result =
                    serde_json::json!({"schema_version":1,"failures":failures,"jobs":results});
                let json = serde_json::to_string_pretty(&result)?;
                if let Some(path) = report {
                    use std::io::Write;
                    let mut file = std::fs::OpenOptions::new()
                        .write(true)
                        .create_new(true)
                        .open(path)?;
                    file.write_all(json.as_bytes())?;
                }
                println!("{json}");
                if failures > 0 {
                    anyhow::bail!("{failures} batch job(s) failed; see results");
                }
            }
            Command::Validate { input } => {
                let result = dnoise::validation::inspect(&input, true)?;
                println!(
                    "Valid type-2 {} acquisition: {} frames, {} points, {} binary bytes",
                    dnoise::detect_acquisition(&input)?,
                    result.frames,
                    result.points,
                    result.binary_bytes
                );
            }
            Command::Metadata { input } => match dnoise::provenance::read(&input)? {
                Some(history) => println!("{}", serde_json::to_string_pretty(&history)?),
                None => println!(
                    "No dnoise history found (older versions or removed sidecars cannot be detected)."
                ),
            },
        }
        Ok(())
    })())
}
