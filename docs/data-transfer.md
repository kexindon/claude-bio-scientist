# Moving the data to the lab server over a cable

101 GB over a Thunderbolt/USB4 cable takes roughly 10–20 minutes. Over Wi-Fi
it takes hours, so use the cable.

Measured on 2026-10-09, source `~/Documents/GitHub/` on the M3:

| Tree | Size | Files |
| --- | --- | --- |
| `Database/` | 72 GB | 1,234 |
| `Output/` | 29 GB | 1,650 |
| `Input/` | 628 MB | 81 |

One file, `Database/Pathogenecity/dbNSFP/dbNSFP5.2a_grch37.gz`, is **40 GB —
40% of everything** — and is referenced by only 2 files in the whole codebase.
Skipping it brings the transfer down to 61 GB. It is re-downloadable from the
dbNSFP site if ever needed. Your call; the commands below include it.

---

## 0. Before you start

On the **M1 Max**, make sure the destination exists and you own it:

```bash
sudo mkdir -p /Users/Shared/lab-data
sudo chown "$(whoami):staff" /Users/Shared/lab-data
df -h /Users/Shared          # confirm >110 GB free
```

On the **M3**, note what you are sending, to compare against later:

```bash
cd ~/Documents/GitHub
for d in Database Output Input; do
  echo "$d  $(find $d -type f | wc -l) files  $(du -sh $d | cut -f1)"
done
```

Write those numbers down.

## 1. Connect the cable and find the direct link

Plug the two Macs together with a **Thunderbolt 4 / USB4 cable**. A
charge-only USB-C cable will not carry data.

On the **M1 Max**:

```bash
ifconfig | grep -B3 "inet 169.254" | grep "^en"
ifconfig | grep "inet 169.254"
```

You want the `169.254.x.x` address — that is the direct cable link. Write it
down; the commands below call it `<CABLE_IP>`.

**Do not use `m1max` for the transfer.** That name routes over Tailscale,
which is an encrypted tunnel and much slower than the cable. Use the raw
`169.254.x.x` address so the bytes go down the wire.

If no `169.254` address appears, check **System Settings → Network** on both
machines for a new Thunderbolt Bridge interface and that it is active.

## 2. Transfer

Run these from the **M3**, one at a time. Start with `Input` — it is small and
proves the path works before you commit to 72 GB.

```bash
cd ~/Documents/GitHub

# smallest first, as a smoke test (628 MB)
rsync -avh --progress --partial \
  Input/ kexin@<CABLE_IP>:/Users/Shared/lab-data/Input/

# 29 GB
rsync -avh --progress --partial \
  Output/ kexin@<CABLE_IP>:/Users/Shared/lab-data/Output/

# 72 GB — the long one
rsync -avh --progress --partial \
  Database/ kexin@<CABLE_IP>:/Users/Shared/lab-data/Database/
```

To skip the 40 GB dbNSFP file, add `--exclude` to the last command:

```bash
rsync -avh --progress --partial --exclude 'Pathogenecity/dbNSFP/' \
  Database/ kexin@<CABLE_IP>:/Users/Shared/lab-data/Database/
```

**The trailing slashes matter.** `Database/` means "the contents of", while
`Database` would nest it as `.../Database/Database/`.

`--partial` keeps incomplete files, so an interrupted run resumes instead of
starting over. If it drops, just run the same command again.

## 3. Verify before deleting anything

Run the second rsync pass — it should transfer nothing:

```bash
rsync -avh --dry-run --delete-after \
  Database/ kexin@<CABLE_IP>:/Users/Shared/lab-data/Database/ | tail -5
```

Then compare counts and sizes on both ends:

```bash
# on the M3
cd ~/Documents/GitHub
for d in Database Output Input; do
  echo "$d  $(find $d -type f | wc -l) files  $(du -sh $d | cut -f1)"
done

# on the M1 Max
cd /Users/Shared/lab-data
for d in Database Output Input; do
  echo "$d  $(find $d -type f | wc -l) files  $(du -sh $d | cut -f1)"
done
```

The file counts must match. Then checksum the files that matter most — the
ones every analysis loads:

```bash
# run on BOTH machines, compare the output
find . -name "*.db" -o -name "*genomic.fna.gz" | sort | xargs shasum -a 256
```

**Do not delete anything from the M3 until those checksums match.** Right now
this is your only copy.

## 4. Make it read-only to other users

On the **M1 Max**, so another account's script cannot corrupt the reference
data everyone depends on:

```bash
chmod -R a+rX,go-w /Users/Shared/lab-data
ls -ld /Users/Shared/lab-data
```

## 5. Point the environment at it

On the **M1 Max**, in `~/.zshrc`:

```bash
export LAB_DATA_ROOT="/Users/Shared/lab-data/Database"
export LAB_OUTPUT_ROOT="/Users/Shared/lab-data/Output"
export LAB_INPUT_ROOT="/Users/Shared/lab-data/Input"
```

Then `source ~/.zshrc` and check a real file resolves:

```bash
ls -l "$LAB_DATA_ROOT/Genecode/release_m33_GRCm39/" | head
python3 -c "import os; print(os.path.exists(os.environ['LAB_DATA_ROOT']))"
```

Every repo's code now reads these variables, so the same notebook runs on both
machines with no edits. See [`data-layout.md`](data-layout.md).

## 6. Afterwards

Keep the M3 copy until the server has been used for real work for a week or
two. The M3 is at 93% disk, so freeing 101 GB is tempting — but a verified
second copy is the only backup this data has. When you do delete it, get an
external drive first so there are still two copies somewhere.
