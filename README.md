# PrimerDesigner

PrimerDesigner is an open-source web workbench for primer design and thermodynamic review. It combines a SvelteKit frontend, a FastAPI backend, Primer3 primer design, primer3-py dimer/Tm analysis, seqfold hairpin prediction, and a custom Mathews/Turner hairpin second opinion.

## Current Product

- Design PCR primer pairs from a template sequence or FASTA paste using exact-bound, required-region, or exploratory internal amplicon modes.
- Analyze pasted primer sets for Tm, GC, hairpins, self-dimers, and cross-dimers.
- Validate Golden Gate overhangs against uniqueness, reverse-complement collisions, palindromes, NEB high-fidelity sets, and internal Type IIS sites.
- Export design and analysis results as CSV.
- Keep advanced reaction conditions available without making the first screen complex.

## Architecture

```text
PrimerDesigner/
  primerdesignr/          # pure Python analysis/design engine
  api/                    # FastAPI routes and web schemas
  frontend/               # SvelteKit web app
  tests/                  # launch workflow regression tests
  Dockerfile              # backend container for hosted API
  railway.json            # Railway backend deployment
```

Hosted deployment is intentionally split:

- Frontend: SvelteKit on Vercel.
- Backend: FastAPI container on Railway, Fly, Render, or similar.

Self-hosting as a single web app can be added later, but the current split keeps the scientific Python runtime off the frontend host and makes the public demo simple.

## Local Development

Backend:

```bash
cd api
pip install -r requirements.txt
uvicorn main:app --reload --port 8000
```

Frontend:

```bash
cd frontend
npm install
npm run dev
```

Create `frontend/.env`:

```bash
VITE_API_URL=http://localhost:8000
```

API docs are available at `http://localhost:8000/docs`.

## Environment

Backend:

```bash
CORS_ORIGINS=https://your-frontend.vercel.app,http://localhost:5173
```

Frontend:

```bash
VITE_API_URL=https://your-api.up.railway.app
```

Production must set both values. Without `VITE_API_URL`, the production frontend will call same-origin API paths, which only works if you add a reverse proxy.

## API

| Method | Endpoint | Purpose |
| --- | --- | --- |
| `GET` | `/health` | Service health check |
| `POST` | `/design` | Design ranked PCR primer candidates from a template |
| `POST` | `/analyze` | Analyze 1-24 existing primers and cross-dimers |
| `POST` | `/pair` | Analyze a forward/reverse pair |
| `POST` | `/golden-gate` | Validate Golden Gate overhangs and scan internal sites |

## Validation

```bash
python3 -m unittest discover
python3 -m compileall primerdesignr api
cd frontend && npm run build
```

## Launch Notes

- The backend is CORS allowlisted. Add your frontend origin in `CORS_ORIGINS`.
- Input limits are enforced server-side: 24 primers per analysis, 200 nt per primer, 20 kb per design template.
- Python and frontend dependency versions are pinned or locked for repeatable deploys.
- `primer3-py` wraps Primer3. Review Primer3/GPL licensing before choosing the repository license and before embedding the backend in non-GPL products.
