# Publishing pyranges0

The GitHub repository is `pyranges/pyranges0`; its PyPI package is `pyranges`.

## One-time PyPI setup

On the [pyranges publishing settings page](https://pypi.org/manage/project/pyranges/settings/publishing/),
add a GitHub Trusted Publisher with these values:

| Field | Value |
|---|---|
| Owner | `pyranges` |
| Repository name | `pyranges0` |
| Workflow name | `build_and_upload_wheels.yml` |
| Environment name | `pypi` |

Use only the workflow filename, without `.github/workflows/`.
No PyPI username, account password, or API-token secret is needed by this workflow.
See the [PyPI setup instructions](https://docs.pypi.org/trusted-publishers/adding-a-publisher/).

## Publish a release

1. Update the version in `pyproject.toml` and the changelog, then merge the changes into `master`.
2. Wait for CI and the distribution build to pass.
3. Open GitHub Actions, select **Build and upload to PyPI**, and choose **Run workflow** on `master`.
4. Verify that both the wheel and source distribution are listed on PyPI.

Pull requests build and validate the distributions without publishing them.
Only a manual run on `master` can execute the publishing job.
The build job has read-only repository access; only the separate publishing job
can obtain the short-lived PyPI publishing credential through GitHub OIDC.

If a release fails before any files are uploaded, fix the workflow and start a new
run from the updated `master`. Rerunning an old failed run uses its old workflow.
If files have already been uploaded, inspect PyPI before retrying: published files
cannot be replaced by rebuilding the same version.
