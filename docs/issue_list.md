## Updated Issue List (Rank-Ordered)

### Critical (Correctness)

1. ~~Z-direction halo bug in integer version~~ — ✓ Fixed

2. ~~Inconsistent array indexing throughout~~ — ✓ Standardized

3. ~~Halo/boundary ordering semantics undefined~~ — ✓ Documented, mutually exclusive by design

### Architectural (Maintainability/Extensibility)

4. ~~Hidden global state~~ — ✓ Fixed, `dirichlet_values` and `bc_types` now passed explicitly

5. **Boundary config not self-contained** — no extensibility path for custom BCs (function pointers, etc.)

6. ~~Code duplication for types~~ — Accepted as technical debt

### Testing Gaps

7. **No integer version of tests** — only `check_halo_real` and `check_boundary_real` exist

8. **Tolerance below SP epsilon** — `1.0e-10_sp` meaningless for single precision

### Performance

9. **No non-blocking communication** — `MPI_Sendrecv` serializes

10. **Stack-allocated buffers every call** — repeated allocation in hot path

### Minor

11. **`domain%abort` signature issue** — `ierr` usage incorrect

12. **DEBUG output garbled** — no rank ordering (your test output demonstrates this nicely)
