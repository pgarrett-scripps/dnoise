// Embeds the Windows icon resource. winresource is only a build-dependency on
// Windows hosts, so all references to it must be compiled out elsewhere. The
// icon is excluded from the crates.io package, so a missing file is a no-op
// and `cargo install dnoise` behaves as before.

#[cfg(windows)]
fn embed_icon() {
    if std::env::var("CARGO_CFG_TARGET_OS").as_deref() != Ok("windows") {
        return;
    }
    let icon = "assets/dnoise.ico";
    if std::path::Path::new(icon).exists() {
        let mut res = winresource::WindowsResource::new();
        res.set_icon(icon);
        if let Err(e) = res.compile() {
            println!("cargo:warning=icon embedding failed: {e}");
        }
    }
}

#[cfg(not(windows))]
fn embed_icon() {}

fn main() {
    embed_icon();
    println!("cargo:rerun-if-changed=assets/dnoise.ico");
}
