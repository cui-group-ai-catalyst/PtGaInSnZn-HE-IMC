# UMA Checkpoint Setup
## Reviewer-Facing Setup Notes for Panels C and G

## 1. Source

The `UMA` checkpoint used by this project comes from:

- Hugging Face model repository: `facebook/UMA`
- URL: `https://huggingface.co/facebook/UMA`

This repository is gated and currently distributed under:

- `FAIR Chemistry License v1`

The package therefore should not imply that the checkpoint is a project-owned
asset. It is an external model dependency that must be obtained under the
upstream license terms.

## 2. Checkpoint used in this project

The project uses:

- checkpoint name: `uma-s-1p1.pt`

The Hugging Face model page lists the MD5 checksum for `uma-s-1p1.pt` as:

- `36a2f071350be0ee4c15e7ebdd16dde1`

After downloading the checkpoint, place it at the default `huggingface_hub` cache
location (or any location of your choice — the panel scripts read `UMA_CHECKPOINT_PATH`
from the environment, falling back to the default cache):

```
# Linux / macOS (default)
~/.cache/huggingface/hub/models--facebook--UMA/snapshots/38529caa2c51a9a8a0d71f0b56b79ac33bc9eceb/checkpoints/uma-s-1p1.pt

# Windows (default)
%USERPROFILE%\.cache\huggingface\hub\models--facebook--UMA\snapshots\38529caa2c51a9a8a0d71f0b56b79ac33bc9eceb\checkpoints\uma-s-1p1.pt

# Or: arbitrary location
export UMA_CHECKPOINT_PATH=/your/local/path/uma-s-1p1.pt
```

After placing the file, verify the MD5 matches the Hugging Face listing
(`36a2f071350be0ee4c15e7ebdd16dde1`) before running any UMA-dependent panel.

## 3. Why this matters

Panels `C` and `G` depend on this checkpoint for numerical reproduction.

Without an explicit checkpoint route, a reviewer may understand the code path
but still be unable to reproduce the same results from a clean machine.

## 4. Recommended setup statement for reviewers

Reviewers should be instructed to:

1. obtain access to `facebook/UMA` on Hugging Face under the FAIR Chemistry
   License v1;
2. download `uma-s-1p1.pt`;
3. verify that its MD5 checksum is
   `36a2f071350be0ee4c15e7ebdd16dde1`;
4. place it in the expected cache location, or adapt the script path
   explicitly if a different local checkpoint path is used.

## 5. Current script behavior

The current project scripts for `FigC` and `FigG` expect the checkpoint to be
present in the local cache path. They do not yet include an automatic download
or cache-bootstrap routine.

This is acceptable for an internal or reviewer package if documented clearly,
but it should be stated explicitly in the release README.

## 6. Practical release note

The safest wording is:

`Panels C and G require the gated UMA checkpoint uma-s-1p1.pt from the Hugging Face repository facebook/UMA. The checkpoint used in this study has MD5 36a2f071350be0ee4c15e7ebdd16dde1 and is expected at the documented local cache path unless a user edits the script configuration.`
