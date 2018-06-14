# Nix development environment
#
# source <bbp-nixpkgs>/sourcethis.sh
#
# nix-build
# nix-shell
#
with import <nixpkgs> {};
{
  neuroc = neuroc.overrideDerivation (oldAtr: rec {
    version = "DEV_ENV";
    src = ./.;
  });
}
