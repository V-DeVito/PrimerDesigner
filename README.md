# primerdesignr-web

Web interface for primerdesignr. SvelteKit frontend + FastAPI backend.

## Architecture

```
┌──────────────────────┐     ┌──────────────────────┐
│  SvelteKit frontend  │────▶│   FastAPI backend     │
│  Vercel (free)       │     │   Railway (free)      │
│                      │     │                       │
│  - Primer input      │     │  - primer3-py (Tm,    │
│  - Results table     │     │    dimers)             │
│  - Cross-dimer matrix│     │  - seqfold (hairpin)  │
│  - CSV export        │     │  - mathews_hairpin    │
│  - IDT order copy    │     │  - Golden Gate        │
└──────────────────────┘     └──────────────────────┘
```

## Local Development

### Backend

```bash
cd api
pip install -r requirements.txt
uvicorn main:app --reload --port 8000
```

API docs: http://localhost:8000/docs

### Frontend

```bash
cd frontend
npm install
npm run dev
```

App: http://localhost:5173

### Environment

Create `frontend/.env`:
```
VITE_API_URL=http://localhost:8000
```

For production, set this to your Railway URL.

## Deploy

### Backend → Railway

1. Push this repo to GitHub
2. Create a new Railway project, connect the repo
3. Set the root directory to `/` (Dockerfile is at project root)
4. Railway auto-detects the Dockerfile and deploys
5. Note your Railway URL: `https://your-app.up.railway.app`

### Frontend → Vercel

1. Import the repo in Vercel
2. Set root directory to `frontend/`
3. Set environment variable: `VITE_API_URL=https://your-app.up.railway.app`
4. Deploy

### DNS (optional)

If using `primerdesignr.com`:
- Point `primerdesignr.com` → Vercel (frontend)
- Point `api.primerdesignr.com` → Railway (backend)
- Update CORS in `api/main.py` to include your domain

## API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/health` | Healthcheck |
| POST | `/analyze` | Analyze 1-24 primers with cross-dimer matrix |
| POST | `/pair` | Analyze a specific primer pair |
| POST | `/golden-gate` | Validate Golden Gate overhangs |

See http://localhost:8000/docs for full schema.

## Project Structure

```
primerdesignr-web/
├── Dockerfile              # Railway deployment
├── api/
│   ├── main.py             # FastAPI endpoints
│   └── requirements.txt
├── primerdesignr/           # Core engine (Python package)
│   ├── __init__.py
│   ├── thermo.py           # Tm, dimers, hairpin analysis
│   ├── mathews_hairpin.py  # AT-closing stem second opinion
│   └── assembly.py         # Gibson, Golden Gate logic
├── frontend/
│   ├── package.json
│   ├── svelte.config.js
│   ├── vite.config.js
│   ├── src/
│   │   ├── app.html
│   │   ├── app.css          # Design tokens, base styles
│   │   ├── lib/
│   │   │   └── api.js       # API client + parser
│   │   └── routes/
│   │       ├── +layout.svelte
│   │       └── +page.svelte  # Main analyzer UI
│   └── .env
└── tests/
    └── test_amye_primers.py
```
