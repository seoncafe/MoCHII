# Final hydrogen and atom-reader checks

The metal selector must not redirect common H I excitation. Actual registration
file selection and the reader loaded the existing 25-level H I atom even with
an invalid metal selector. The diagnostic returns after the production reader
before numerical table construction. No H I table or scientific product was
regenerated for the path correction. The complete production build passed.

Five installed endpoint/Fe II atoms loaded through the current production
parser. Six malformed upper-endpoint atoms were rejected: negative, NaN,
below the final node, extra values, duplicate, and no collision rows. The
empty atom is rejected by the existing invalid NUPS check before its endpoint
record. Both parser stages prohibit an endpoint model with no collision rows.

The preceding attempt correctly rejected a negative ceiling, but its wrapper
incorrectly expected a NLEVEL prefix; the production reader prints ERROR.
Its failed log and diagnostic executable were preserved. The continuation
checks actual reasons. Its final empty-atom assertion originally expected
the later endpoint error; the earlier NUPS error correctly prevented reading
the invalid model. All six saved failures were checked without rerunning them.
