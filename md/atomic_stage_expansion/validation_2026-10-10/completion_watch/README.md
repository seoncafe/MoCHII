# Calculation completion detection and automatic continuation

The user requested automatic detection so finished calculations are not overlooked.

The reusable calculation controller validates each expected JSON result and runs the next declared command. The event watcher additionally detects completed steps and controller failures, rejects duplicate events, and waits while the parent UI session is active before launching a serialized model execution.

Same-session resume was tested and failed with `thread/resume failed: list_turns is not supported yet (code -32601)`. The supported replacement is a **new execution session** with the documented [handoff](handoff.md); it is not the original UI conversation. Delivery of its output into the original UI was not established.

The replacement was tested end to end: a real numerical check calculated a sum of 55, its controller validated the output, the watcher detected the completion, and a new model session returned `WATCHER_COMPLETION_OK` with exit status zero and a `turn.completed` event. Parent test session: `01a123a5-9acb-7b72-841e-d154778f98d2`. Model execution session: `01a123a6-7496-7de3-8095-e24228df6f7c`. The test watcher was stopped after certification. Evidence is preserved under `smoke/`.

The production watcher runs in tmux session `mochii_completion_watch`, with a 24-hour service limit, a 6 GiB virtual-memory limit, and a one-hour limit for each model execution. Its explicit parent session is `01a11df3-46ae-7a40-a52f-b0bf65878b1e`. Existing reviewed events are baselined; the four completed capacity-study events remain queued for the first continuation. Later registered calculation-controller completions and failures are detected through the configured path patterns.

Inspect [watcher/state.json](watcher/state.json) for actual status. Timestamp directories retain prompts, JSON model events, errors, and `last_response.md`. The model updates `md/atomic_stage_expansion/automatic_continuation_status.md`. A `STOP` file ends watching. Overall completion or exhaustion of independent work pending a required public-design decision also stops the watcher; calculation completion alone is not scientific acceptance.

No external email or chat message is sent. No guarantee of an original-UI notification is made. A stopped, failed, or timed-out watcher requires explicit inspection; uncertain model work is not silently restarted.

## Failure found during the first actual continuation

Completion detection and a model text response passed the initial smoke test,
but that test did not exercise shell tools. The first scientific handoff
failed every execution call with `code-mode host closed its stdout`, despite
a zero model exit and a completed model turn. It therefore performed no
scientific continuation. A separate real-shell probe with the host disabled
also failed with `code-mode host is disabled`; it did not return the requested
`WATCHER_TOOL_OK`. Neither probe certifies executable model continuation.

The main session has placed `STOP` to stop the ineffective watcher. No running
scientific calculation was stopped. The controller capacity studies are
complete, but scientific work remains incomplete and the main session must
continue it directly until a working automatic execution route is verified.
