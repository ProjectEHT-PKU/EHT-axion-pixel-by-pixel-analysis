"""Execute the inference research stage using the provided data and configuration."""

from _stage import execute


def main(argv=None):
    return execute("inference", argv)


if __name__ == "__main__":
    main()
