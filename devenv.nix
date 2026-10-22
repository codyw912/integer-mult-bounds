{ pkgs, ... }:

{
  # Registered research-package replays use Python 3.12 and the standard library.
  packages = [ pkgs.git pkgs.uv ];
  languages.python = {
    enable = true;
    package = pkgs.python312;
  };
}
