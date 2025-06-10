import csv
import subprocess
import os
import sys
import shutil
from collections import Counter

def main():
    base_dir        = os.path.abspath(os.path.dirname(__file__))
    dir_653         = os.path.join(base_dir, "653")
    csv_path        = os.path.join(dir_653, "Table_S1_from_653.csv")
    working_dir     = os.path.join(dir_653, "Working")
    half_dir        = os.path.join(dir_653, "Half_working")
    fake_dir        = os.path.join(dir_653, "Fake_working")
    large_dir       = os.path.join(dir_653, "Working_Large")
    non_working_path= os.path.join(dir_653, "non-working.txt")
    summary_path    = os.path.join(dir_653, "summary.txt")
    errors_path     = os.path.join(dir_653, "errors.txt")

    # ensure folders exist
    for d in (working_dir, half_dir, fake_dir, large_dir):
        os.makedirs(d, exist_ok=True)

    working_count       = 0
    half_working_count  = 0
    fake_working_count  = 0
    large_working_count = 0
    non_working         = []
    error_counter       = Counter()

    with open(csv_path, newline='') as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            if len(row) < 3:
                print(f"Skipping row with too few columns: {row}", flush=True)
                continue

            comp_id, smiles = row[2].strip(), row[1].strip()
            print(f"Attempt ID {comp_id}...", flush=True)

            try:
                # 1) Run main_batch.py inside Working/
                result = subprocess.run(
                    [sys.executable, os.path.join(base_dir, "main_batch.py")],
                    input=f"{comp_id}\n{smiles}\n",
                    text=True,
                    cwd=working_dir,
                    capture_output=True,
                )
                if result.returncode != 0:
                    msg = (result.stdout or result.stderr).strip()
                    raise RuntimeError(msg)

                # 2) Create per‐compound folder
                comp_folder = os.path.join(working_dir, comp_id)
                if os.path.isdir(comp_folder):
                    shutil.rmtree(comp_folder)
                os.makedirs(comp_folder)

                # 3) Move all three outputs: .gro, .itp, .txt
                for ext in ("gro", "itp", "txt"):
                    fname = f"{comp_id}.{ext}"
                    src   = os.path.join(working_dir, fname)
                    if not os.path.exists(src):
                        raise FileNotFoundError(f"Missing expected output: {fname}")
                    shutil.move(src, os.path.join(comp_folder, fname))

                # 4) Read .txt to classify fake (instead of using the .gro file)
                txt_out_path = os.path.join(comp_folder, f"{comp_id}.txt")
                with open(txt_out_path) as tf:
                    txt_lines = tf.readlines()
    
                # “fake” if any line (after the first 6) has more than 10 commas
                is_fake = any(line.count(',') > 8 for line in txt_lines[6:])
    
                # 5) Read .itp to classify half‑working
                itp_path = os.path.join(comp_folder, f"{comp_id}.itp")
                with open(itp_path) as itf:
                    itp_lines = itf.readlines()
                is_half = (not is_fake) and any('z' in line for line in itp_lines[6:])
                
                # still use .gro to detect “large”
                gro_path = os.path.join(comp_folder, f"{comp_id}.gro")
                with open(gro_path) as gf:
                    gro_lines = gf.readlines()
                is_large = (not is_fake) and (not is_half) and len(gro_lines) >= 25

                if is_fake:
                    dest = os.path.join(fake_dir, comp_id)
                    if os.path.isdir(dest): shutil.rmtree(dest)
                    shutil.move(comp_folder, dest)
                    fake_working_count += 1
                    print(f"  • {comp_id}: FAKE → Fake_working/{comp_id}/", flush=True)
                    continue

                if is_half:
                    dest = os.path.join(half_dir, comp_id)
                    if os.path.isdir(dest): shutil.rmtree(dest)
                    shutil.move(comp_folder, dest)
                    half_working_count += 1
                    print(f"  • {comp_id}: HALF → Half_working/{comp_id}/", flush=True)
                    continue

                if is_large:
                    dest = os.path.join(large_dir, comp_id)
                    if os.path.isdir(dest): shutil.rmtree(dest)
                    shutil.copytree(comp_folder, dest)
                    large_working_count += 1
                    working_count      += 1
                    print(f"  • {comp_id}: LARGE → Working_Large/{comp_id}/", flush=True)
                    continue

                # fully working
                working_count += 1
                print(f"  • {comp_id}: WORKING → Working/{comp_id}/", flush=True)

            except Exception as e:
                print(f"  ✗ {comp_id}: ERROR – {e}", flush=True)
                error_counter[str(e)] += 1
                non_working.append(comp_id)
                # clean up any partial folder
                bad = os.path.join(working_dir, comp_id)
                if os.path.isdir(bad):
                    shutil.rmtree(bad)

    # write non‐working list
    with open(non_working_path, "w") as fw:
        for cid in non_working:
            fw.write(f"{cid}\n")

    # write summary
    with open(summary_path, "w") as fs:
        fs.write(f"working count: {working_count}\n")
        fs.write(f"half-working count: {half_working_count}\n")
        fs.write(f"fake-working count: {fake_working_count}\n")
        fs.write(f"large-working count: {large_working_count}\n")
        fs.write(f"non-working count: {len(non_working)}\n")

    # write error summary
    with open(errors_path, "w") as fe:
        fe.write("Error summary (counts per exact message):\n")
        for err, cnt in error_counter.most_common():
            fe.write(f"{cnt} times Error: {err}\n")

    # final report
    print("\nDone.")
    print(f"  • {working_count} fully successful runs   → Working/<ID>/")
    print(f"  • {half_working_count} half‑working runs   → Half_working/<ID>/")
    print(f"  • {fake_working_count} fake‑working runs   → Fake_working/<ID>/")
    print(f"  • {large_working_count} large‑working runs  → Working_Large/<ID>/")
    print(f"  • {len(non_working)} failures             → non-working.txt")
    print("  • summary.txt written.")
    print("  • errors.txt written.")

if __name__ == "__main__":
    main()
