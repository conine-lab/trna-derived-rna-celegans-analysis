# Standard Operating Procedure (SOP): Creating and maintaining analysis repositories

> **What this is:**  
> This document is a **Standard Operating Procedure (SOP)** — a written guide
> describing the agreed-upon, correct way the lab creates, edits, and publishes
> analysis repositories.
>
> This document is **internal lab guidance** and is not part of the published
> methods or analysis pipeline.

This SOP describes the **standard lab process** for creating, documenting,
and maintaining analysis repositories associated with publications.

Follow this SOP unless there is a clear reason not to.

---

## Background and reference material

The practices in this SOP are informed by widely used community standards for
reproducible computational research and Git-based workflows.

Helpful background (optional viewing):

- **How to Make Your Research Code Open Source (YouTube)**  
  https://www.youtube.com/watch?v=r8jQ9hVA2qs

---

## When to create a new repository

Create a new repository when:
- a paper includes custom analysis code
- figures depend on scripted data processing
- code will be cited or archived (e.g. via Zenodo)

Do **not** reuse old repositories across papers.

---

## Repository naming

Repository names should be:
- topic-based (not author-based)
- lowercase
- hyphen-separated

**Good examples**
- `trna-derived-rna-celegans-analysis`
- `small-rna-inheritance-analysis`

**Avoid**
- `galambos-et-al`
- `paper1`
- `analysis`

---

## Creating a new repository (GitHub UI)

1. Create a new **public** repository under the lab organization
2. Add a **LICENSE** (MIT recommended)
3. Optionally add the default README (can be replaced later)

Do not add data or code yet.

---

## Standard repository structure

Every analysis repository must use this structure:

```
repo/
├── scripts/     # analysis code only
├── docs/        # documentation and SOPs
├── examples/    # representative outputs (optional)
├── README.md    # short landing page
├── LICENSE
└── requirements.txt (if applicable)
```

Do not deviate without discussion.

---

## Adding code to the repository

### Before committing scripts

- Clean up filenames (no `_final`, `_v2`, `_updated`)
- Add `-h` help text to all scripts
- Make inputs and outputs explicit
- Remove debugging prints

Scripts should be understandable without opening the paper.

---

### Adding scripts

1. Place scripts in `scripts/`
2. Run each script locally with `-h`
3. Confirm no errors on import

Only finalized scripts should be committed.

---

## Documentation requirements

### Required documentation

Every repository **must include**:

- `README.md` (short, high-level)
- `docs/Code_Environment_Reproducibility.md` (detailed)

Documentation should explain:
- environment setup
- script order
- expected inputs and outputs

Avoid duplicating documentation across files.

---

## Example outputs

If example outputs are included:
- place them in `examples/`
- add `examples/README.md`
- do not treat examples as inputs

Replace example outputs only intentionally and document the change in commit messages.

---

## Version control practices

- Commit early, commit clearly
- One logical change per commit
- Descriptive commit messages required

**Good**
```
Add finalized analysis scripts and documentation
```

**Bad**
```
updates
```

Do not force-push to `main`.

---

## Editing an existing repository

When modifying a repository after initial release:

1. Make changes locally
2. Review `git status` carefully
3. Ensure no OS junk files (e.g. `.DS_Store`)
4. Commit only intentional changes
5. Push to `main`

If changes affect figures or results, update documentation accordingly.

---

## Handling Markdown documentation

- Long documentation belongs in `docs/`
- GitHub only auto-renders `README.md` by default
- Use a short root README to link to longer documentation

Do not rely on users finding documentation by browsing directories.

---

## Preparing for publication or archiving

Before submission:
- repository must be public
- documentation must be complete
- scripts must run without modification

After acceptance:
- create a release tag (e.g. `v1.0.0`)
- link GitHub to Zenodo
- add DOI to `CITATION.cff`

---

## Guiding principles

- Clarity > cleverness
- Consistency > convenience
- If something is confusing, document it
- If it cannot be documented, it is not ready

This SOP exists to make repository creation boring, reproducible, and durable.
