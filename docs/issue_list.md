## Issue List (Rank-Ordered)

### Critical (Correctness)

1. **Z-direction halo bug in integer version** — sends/receives wrong indices - Fixed

2. **Inconsistent array indexing throughout** — mixing `lbound`, hardcoded `1`, and passed sizes with 0-based arrays (`0:nx+1`). Boundary module writes to index `1`, but that's interior, not the ghost cell at index `0`. Tests may pass for wrong reasons. - Fixed

3. **Halo/boundary ordering semantics undefined**: Neumann BC copies from cell 2→1, but halo exchange also writes to cell - Fixed
  - which overwrites which? What's the intended cell layout?

### Architectural (Maintainability/Extensibility)

4. **Hidden global state** — `dirichlet_value` and `boundaries` pulled from module, not passed explicitly. Tests also depend on global `boundaries`.

5. **Boundary config not self-contained** — no per-face Dirichlet values, no extensibility path

6. **Code duplication for types** — **Accepted as technical debt** with strong test coverage

### Testing Gaps

7. **No integer version of tests** — only `check_halo_real` and `check_boundary_real` exist

8. **Tolerance below SP epsilon** — `1.0e-10_sp` is meaningless for single precision; effectively testing exact equality

### Performance

9. **No non-blocking communication** — `MPI_Sendrecv` serializes; overlap potential lost

10. **Stack-allocated buffers every call** — repeated allocation overhead in hot path

### Minor

11. **`domain%abort` signature issue** — `ierr` usage incorrect

12. **DEBUG output garbled** — no rank ordering for multi-process output
