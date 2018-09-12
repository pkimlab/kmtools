from pathlib import Path


def package_file_refs(obj) -> None:
    for key, value in vars(obj).items():
        if key.endswith("_dir") and isinstance(value, Path):
            setattr(obj, key, value.as_posix() if key != "temp_dir" else "")
        if key.endswith("_file") and isinstance(value, Path):
            with value.open("rt") as fin:
                setattr(obj, key, value.name + "\n" + fin.read())


def unpackage_file_refs(obj, temp_dir) -> None:
    for key, value in vars(obj).items():
        if key.endswith("_dir") and isinstance(value, str):
            setattr(obj, key, Path(value) if key != "temp_dir" else temp_dir)
        if key.endswith("_file") and isinstance(value, str):
            name, _, data = value.partition("\n")
            temp_file = temp_dir.joinpath(name)
            with temp_file.open("wt") as fout:
                fout.write(data)
            setattr(obj, key, temp_file)
