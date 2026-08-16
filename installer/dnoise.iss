; dnoise Windows installer (Inno Setup 6).
;
; Built by .github/workflows/release.yml on windows-latest via:
;   ISCC.exe /DAppVersion=<version> /Odist installer\dnoise.iss
;
; Per-user install (no admin prompt): files under %LOCALAPPDATA%\Programs\dnoise,
; registry writes under HKCU only. The optional "contextmenu" task adds a
; "Denoise with dnoise" entry to the Explorer right-click menu of .d folders
; (on Windows 11 it appears under "Show more options"). The optional
; "addtopath" task appends the install dir to the per-user PATH.

#ifndef AppVersion
  #define AppVersion "0.0.0-dev"
#endif

[Setup]
AppId={{3C9C3286-BE97-4BD4-9BED-8FBDA932212E}
AppName=dnoise
AppVersion={#AppVersion}
AppPublisher=Patrick Garrett
AppPublisherURL=https://github.com/pgarrett-scripps/dnoise
AppSupportURL=https://github.com/pgarrett-scripps/dnoise/issues
PrivilegesRequired=lowest
DefaultDirName={autopf}\dnoise
DisableProgramGroupPage=yes
OutputBaseFilename=dnoise-setup-{#AppVersion}-windows-x86_64
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
Compression=lzma2
SolidCompression=yes
SetupIconFile=..\assets\dnoise.ico
UninstallDisplayIcon={app}\dnoise-gui.exe
ChangesEnvironment=yes
LicenseFile=..\LICENSE

[Tasks]
Name: "contextmenu"; Description: "Add ""Denoise with dnoise"" to the right-click menu of .d folders"
Name: "addtopath"; Description: "Add dnoise to PATH (command-line use)"; Flags: unchecked

[Files]
Source: "..\target\release\dnoise.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\target\release\dnoise-gui.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\README.md"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\LICENSE"; DestDir: "{app}"; Flags: ignoreversion

[Icons]
Name: "{autoprograms}\dnoise"; Filename: "{app}\dnoise-gui.exe"

[Registry]
; Explorer context menu on .d folders. AppliesTo (an AQS query) hides the entry
; on other folders; if a Windows build ignores it, a wrong click is harmless
; because dnoise-gui rejects non-.d folders with a log note.
Root: HKCU; Subkey: "Software\Classes\Directory\shell\dnoise"; \
  ValueType: string; ValueData: "Denoise with dnoise"; \
  Tasks: contextmenu; Flags: uninsdeletekey
Root: HKCU; Subkey: "Software\Classes\Directory\shell\dnoise"; \
  ValueType: string; ValueName: "AppliesTo"; ValueData: "System.FileName:""*.d"""; \
  Tasks: contextmenu
Root: HKCU; Subkey: "Software\Classes\Directory\shell\dnoise"; \
  ValueType: string; ValueName: "Icon"; ValueData: """{app}\dnoise-gui.exe"""; \
  Tasks: contextmenu
Root: HKCU; Subkey: "Software\Classes\Directory\shell\dnoise\command"; \
  ValueType: string; ValueData: """{app}\dnoise-gui.exe"" ""%1"""; \
  Tasks: contextmenu
; Per-user PATH append, only if not already present.
Root: HKCU; Subkey: "Environment"; ValueType: expandsz; ValueName: "Path"; \
  ValueData: "{olddata};{app}"; Tasks: addtopath; \
  Check: NeedsAddPath(ExpandConstant('{app}'))

[Code]
// True when {app} is not already on the per-user PATH.
function NeedsAddPath(Param: string): Boolean;
var
  Path: string;
begin
  if not RegQueryStringValue(HKCU, 'Environment', 'Path', Path) then
  begin
    Result := True;
    exit;
  end;
  Result := Pos(';' + Uppercase(Param) + ';', ';' + Uppercase(Path) + ';') = 0;
end;

// Strip {app} from the per-user PATH on uninstall.
procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
var
  Path, App: string;
  P: Integer;
begin
  if CurUninstallStep <> usPostUninstall then exit;
  if not RegQueryStringValue(HKCU, 'Environment', 'Path', Path) then exit;
  App := ExpandConstant('{app}');
  P := Pos(';' + Uppercase(App), ';' + Uppercase(Path));
  if P > 0 then
  begin
    Delete(Path, P, Length(App) + 1);
    RegWriteExpandStringValue(HKCU, 'Environment', 'Path', Path);
  end;
end;
