# Galaxy wrapper deployment

The wrapper is published to the Galaxy Tool Shed automatically. 

## Versioning

The tool version is set in `macros.xml`:

```
<TOOL_VERSION>+galaxy<VERSION_SUFFIX>     e.g. 1.2.0+galaxy1
```

- `@TOOL_VERSION@`: the BeeGees version on Bioconda.
- `@VERSION_SUFFIX@`: the wrapper revision for that BeeGees version. It
  resets to `0` whenever `@TOOL_VERSION@` changes.

## Workflows

| Workflow | Trigger | What it does |
|---|---|---|
| `galaxy_autoupdate.yml` | Daily, or manual | Opens a PR bumping `@TOOL_VERSION@` when Bioconda has a newer BeeGees |
| `galaxy_tag.yml` | `macros.xml` changes on `main` | Pushes tag `galaxy-<TOOL_VERSION>+galaxy<VERSION_SUFFIX>` if it doesn't exist yet |
| `toolshed_deploy.yml` | Any `galaxy-*` tag | Lints, tests and uploads to the Tool Shed |

## Releasing

### Wrapper-only change

1. Edit the wrapper files.
2. Bump `@VERSION_SUFFIX@` in `macros.xml`.
3. Merge to `main`. The tag and Tool Shed upload follow automatically.

Changes merged without a suffix bump are not released.

### New BeeGees release on Bioconda

1. The daily job opens a PR on `galaxy-autoupdate`.
2. If the new release needs wrapper changes, commit them to
   `galaxy-autoupdate` and merge before the next daily run (06:00 UTC).
   Until the PR is merged, that run may overwrite the branch and drop your
   commits. For fixes that take longer, branch off `galaxy-autoupdate`, open
   your own PR, and close the bot's PR once yours is merged.
3. Merge to `main`. The tag and Tool Shed upload follow automatically.

## Requirements

The repository secret `GALAXY_AUTOUPDATE_TOKEN` (a PAT) must be set. Without
it, the autoupdate PR won't run CI and the tag won't trigger the deploy.
