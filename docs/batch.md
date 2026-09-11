# Batch workflows

Save a queue with the GUI's **Save batch…** button, or create a JSON manifest:

```json
{
  "schema_version": 1,
  "jobs": [
    {"input": "raw/A.d", "output": "processed/A.d", "config": {}},
    {"input": "raw/B.d", "output": "processed/B.d",
     "config": {"halo": true, "frame_batch_size": 256}}
  ]
}
```

Paths are relative to the manifest, regardless of the working directory. Create
the `processed` directory first. Each job's `config` uses the TOML field names;
omitted values use the usual CLI defaults.

```sh
dnoise batch jobs.json --report results.json
dnoise batch jobs.json --dry-run
```

Writing runs preflight the whole queue for invalid inputs, settings, and
overlapping destinations before processing starts. Existing output requires
`--force` and is replaced only after that job succeeds. Runtime failures are
recorded per file while later jobs continue. The command exits nonzero if any
job fails; stdout and the optional report contain structured results.

For a retry, create a manifest containing only failed jobs. The GUI's
**Retry unfinished** does this for the current session. Each output also stores
its own completion history and reusable configuration.

For a scheduler, run one manifest per job allocation and keep stdout/stderr as
job logs. Avoid putting the same destination in concurrent allocations. Python
bindings and additional input formats remain separate future extensions.
