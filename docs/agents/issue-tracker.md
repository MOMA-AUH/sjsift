# Issue tracker

Publish specs and tickets as GitHub Issues in MOMA-AUH/sjsift using gh with --repo MOMA-AUH/sjsift.

Read an issue's full body, labels, and comments before acting on it. List existing issues before publishing to avoid duplicates. Use an exact temporary body file with --body-file for multiline issue bodies and comments.

Apply the triage roles defined in triage-labels.md. Ready-for-agent describes specification readiness; respect explicit user instructions about when implementation may begin.

For ticket breakdowns, reference the parent spec in each ticket. Record blockers in the ticket body and use GitHub's native blocking relationships when available. Add a dependency through POST repos/MOMA-AUH/sjsift/issues/<blocked-number>/dependencies/blocked_by, passing the blocker's numeric database ID as issue_id; retrieve that ID from GET repos/MOMA-AUH/sjsift/issues/<blocker-number>. Use textual blocker references if native relationships are unavailable. Preserve the parent issue when publishing tickets.

PRs as a request surface: no.
