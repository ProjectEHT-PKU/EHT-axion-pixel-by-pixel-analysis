"""Execute the raptor research stage using the provided data and configuration."""

from _stage import execute


def main(argv=None):
    return execute("raptor", argv)


if __name__ == "__main__":
    main()
