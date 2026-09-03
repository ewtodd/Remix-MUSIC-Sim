{
  description = "Remix-MUSIC-Sim: catima-based, multi-threaded fork of the ANL MUSIC detector simulator";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/ced43465ad23b2fdea055be721e79895cbf96c28";
    flake-utils.url = "github:numtide/flake-utils";
    catima-src = {
      url = "github:hrosiak/catima/75d22b260ed921f2e6d1c257ca82cf25bcd9f906";
      flake = false;
    };
    # SRIM under wine, for the stopping tables that stopping = "srim" and
    # "mean" read. Its nixpkgs is left alone: the wine wrapper is pinned to
    # what that flake tested.
    srim-nix.url = "github:ewtodd/SRIM-nix";
  };
  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      catima-src,
      srim-nix,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        version = "26.9.2";
        make-srim-table = srim-nix.packages.${system}.make-srim-table;
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
            tomlplusplus
          ];
          buildPhase = ''
            make CATIMA_PREFIX=${catima} VERSION=${version} ASSETS_DIR_OUT=$out/assets \
              SRIM_TABLE_BIN=${make-srim-table}/bin/make-srim-table
          '';
          installPhase = ''
            mkdir -p $out/bin $out/assets
            cp musicsim $out/bin/
            # srim-cache generates the SRIM tables a control file needs by
            # driving SRIM under wine through SRIM-nix's make-srim-table, whose
            # store path it carries; still a deliberate, occasional step, not
            # part of a simulation run.
            cp srim-cache $out/bin/
            cp -r assets/. $out/assets/
          '';
        };
      in
      {
        packages.catima = catima;
        packages.default = musicsim;
        devShells.default = pkgs.mkShell {
          inputsFrom = [ musicsim ];
          nativeBuildInputs = with pkgs; [
            clang-tools
            taplo
            make-srim-table
          ];
          shellHook = ''
            export CATIMA_PREFIX=${catima}
            export MUSICSIM_VERSION=${version}
            export SRIM_TABLE_BIN=${make-srim-table}/bin/make-srim-table
            export CPLUS_INCLUDE_PATH="$PWD/include''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
            export ROOT_INCLUDE_PATH="$PWD/include:${pkgs.root}/include"
            export LD_LIBRARY_PATH="$PWD/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
            # Use .githooks/ for pre-commit formatting (clang-format + taplo).
            git config --local core.hooksPath .githooks 2>/dev/null || true
          '';
        };
      }
    );
}
