fn main() {
    for path in ["src", "Cargo.toml", "Cargo.lock", ".git/HEAD", ".git/index"] {
        println!("cargo:rerun-if-changed={path}");
    }
    let revision = std::process::Command::new("git")
        .args(["describe", "--always", "--dirty"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_owned())
        .unwrap_or_else(|| "unavailable (source archive)".into());
    println!("cargo:rustc-env=DNOISE_BUILD_REVISION={revision}");
}
