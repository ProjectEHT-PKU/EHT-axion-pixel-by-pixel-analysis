"""Build a pinned RAPTOR through its official setup.sh/makefile, never ray-trace."""

from pathlib import Path
import sys

sys.path.insert(0, str(next(p for p in Path(__file__).resolve().parents if (p / "eht_axion").is_dir())))

from eht_axion.raptor.options import options, save_request
from eht_axion.raptor.build import build_raptor


def main(argv=None):
    args = options("build", argv)
    save_request(args)
    output = build_raptor(args.source_root, args.run_root / "build",
                          args.config_data["raptor"]["build"], args.config_data["raptor"]["source"])
    print(f"Verified build: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
