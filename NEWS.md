<!--
This file should document all significant changes brought about by each new
release.

All changes related to a particular release should be collected under a heading
specifying the version number of that release, such as
"# Changes in BioCroValidation Version 2.0.0". The individual changes should be
listed as bullet points and categorized under "## Major Changes",
"## Minor Changes", or "## Bug Fixes" following the major.minor.patch structure
of semantic versioning, or variants of these such as
"## Minor User-Facing Changes".

To facilitate this, when a feature on a feature branch is completed and a pull
request is being prepared, a new section should be added at the top of this file
under the heading "# UNRELEASED"; it should list all the important changes made
on the feature branch.

Then, when it comes time to merge the feature branch into `develop`, the new
"# UNRELEASED" section is transferred into the `develop` branch's version of
NEWS.md, or, if the `develop` branch already has an "# UNRELEASED" section in
its version of NEWS.md, the feature branch's "# UNRELEASED" section will be
integrated into the one on the `develop` branch. (This process of integrating
the two "# UNRELEASED" sections will likely be part of resolving an inevitable
merge conflict.)

Finally, when a new release is made, "# UNRELEASED" should be replaced by a
heading with the new version number, such as
"# Changes in BioCroValidation Version 2.0.0". This section will combine the
draft release notes for all features that have been added since the previous
release.

In the case of a hotfix, a short section headed by the new release number should
be directly added to this file to describe the related changes.
-->

# UNRELEASED

# Changes in BioCroValidation Version 0.2.0 (2025-05-23)

- Added 2002 and 2005 SoyFACE biomass and standard deviation data.

- Added several new functions: `objective_function`, `update_model`, and
  `bounds_table`.

- Added two new vignettes: a "Getting Started" article (`BioCroValidation.Rmd`)
  and a user guide illustrating how to perform a model parameterization
  (`parameterizing_soybean_biocro.Rmd`).

# Changes in BioCroValidation Version 0.1.0

- This is the first version of BioCroValidation. At this point, the package is
  in a state of rapid development, and not all changes will be described here.

- We are reserving version `1.0.0` for a more stable and complete future
  release; until then, major changes should only increase the minor version
  number.
