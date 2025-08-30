#!/usr/bin/env python3
import json
import os
import shlex

def main() -> None:
    sources_str = os.environ.get("SOURCES", "").strip()
    cc = os.environ.get("CC", "clang++").strip()
    cflags = os.environ.get("CFLAGS", "").strip()
    includes = os.environ.get("INCLUDES", "").strip()
    obj_dir = os.environ.get("OBJ_DIR", "obj").strip()

    sources = [s for s in shlex.split(sources_str) if s.endswith('.cpp')]
    directory = os.getcwd()
    entries = []
    for src in sources:
        obj = os.path.join(obj_dir, src.replace('.cpp', '.o'))
        cmd = f"{cc} {cflags} {includes} -c {src} -o {obj}"
        entries.append({
            "directory": directory,
            "file": src,
            "output": obj,
            "command": cmd,
        })

    os.makedirs(os.path.dirname('compile_commands.json') or '.', exist_ok=True)
    with open('compile_commands.json', 'w') as f:
        json.dump(entries, f, indent=2)
    print(f"Wrote compile_commands.json with {len(entries)} entries")

if __name__ == "__main__":
    main()



