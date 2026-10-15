"""Build a pinned RAPTOR through its official setup.sh/makefile, never ray-trace."""
import _bootstrap  # noqa: F401
from _producer import options, save_request
from eht_axion.producer.build import build_raptor


def main(argv=None):
    args = options("build", argv)
    save_request(args)
    output = build_raptor(args.source_root, args.run_root / "build",
                          args.config_data["producer"]["build"], args.expected_commit)
    print(f"Verified build: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
