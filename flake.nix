{
  description = "Remix-MUSIC-Sim: catima-based, multi-threaded fork of the ANL MUSIC detector simulator";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    catima-src = {
      url = "github:hrosiak/catima/75d22b260ed921f2e6d1c257ca82cf25bcd9f906";
      flake = false;
    };
  };
  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      catima-src,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        version = "26.5.12";
        catima = pkgs.stdenv.mkDerivation {
          pname = "catima";
          version = "75d22b2";
          src = catima-src;
          nativeBuildInputs = with pkgs; [ cmake ];
          cmakeFlags = [ "-DBUILD_SHARED_LIBS=ON" ];
        };
        musicsim = pkgs.stdenv.mkDerivation {
          pname = "remix-music-sim";
          inherit version;
          src = ./.;
          nativeBuildInputs = with pkgs; [
            pkg-config
            gnumake
          ];
          buildInputs = with pkgs; [
            root
            catima
          ];
          buildPhase = ''
            make CATIMA_PREFIX=${catima} VERSION=${version}
          '';
          installPhase = ''
            mkdir -p $out/bin
            cp musicsim $out/bin/
          '';
        };
      in
      {
        packages.catima = catima;
        packages.default = musicsim;
        devShells.default = pkgs.mkShell {
          inputsFrom = [ musicsim ];
          nativeBuildInputs = with pkgs; [ clang-tools ];
          shellHook = ''
            export CATIMA_PREFIX=${catima}
            export MUSICSIM_VERSION=${version}
          '';
        };
      }
    );
}
