
# 🧬 Mutalyze

AI-Powered Gene Mutation Impact Predictor

## Overview

Mutalyze is a cloud-based tool that predicts the impact of gene mutations using AI models and public bioinformatics databases. It helps researchers, students, and biotech enthusiasts analyze genetic variations easily, providing confidence scores, protein structure visualization, and prediction history.

## Features

- 🔍 Predict mutation impact (e.g., "Likely Pathogenic")
- 📊 Display confidence scores and external references
- 🗄️ Store prediction history in MongoDB
- 🧬 Visualize protein structure (AlphaFold)
- 📈 Batch analysis via CSV upload
- 🚀 Fast and lightweight deployment
- 📱 User-friendly Streamlit interface

## Architecture

- **Frontend**: Streamlit (Python-based interactive UI)
- **Backend**: FastAPI (REST API for prediction and data)
- **Database**: MongoDB Atlas (for prediction history)
- **ML Model**: RandomForest (scikit-learn, joblib)
- **External APIs**: UniProt, AlphaFold

The frontend communicates with the backend via REST API calls. The backend handles prediction logic, model inference, and data retrieval from external sources.

## Getting Started

### 1. Clone the repository
```bash
git clone https://github.com/yourusername/mutalyze.git
cd mutalyze
```

### 2. Install dependencies
```bash
pip install -r requirements.txt
```

### 3. Set up environment variables
Create a `.env` file in the root directory with:
```
MONGODB_URI=your_mongodb_uri
```

- `MONGODB_URI` is required for storing prediction history. If not set, history will not be persisted.

### 4. Run the application
To start both backend (FastAPI) and frontend (Streamlit), run:
```bash
python run.py
```
- Backend: http://localhost:8000
- Frontend: http://localhost:8501

Alternatively, you can run only the frontend (for demo/testing):
```bash
streamlit run app.py
```

## Usage

1. **Login** with demo credentials (demo@mutalyze.com / demo123)
2. Enter gene name (e.g., BRCA1)
3. Input mutation (e.g., A123T)
4. Click "Predict Impact" to view results
5. Optionally, upload a CSV for batch analysis
6. View prediction history and protein structure visualization

### Example Input
- Gene: `BRCA1`
- Mutation: `A123T`

### Example Output
- Predicted Impact: Likely Pathogenic
- Confidence Score: 0.92
- Protein Name: BRCA1 protein
- Protein Domain: BRCT domain
- Conservation Score: 0.85
- Clinical Significance: Pathogenic
- UniProt ID: P38398
- Structure Visualization: [AlphaFold link]

## API Endpoints

- `POST /predict` — Predict impact for a single mutation
- `POST /batch-predict` — Predict impact for a batch of mutations
- `GET /history` — Retrieve prediction history
- `GET /protein/{uniprot_id}` — Get protein info by UniProt ID

## File Structure

```
MUTALYZE/
  ├── api/                # FastAPI backend
  │   ├── main.py         # API entry point
  │   └── utils/
  │       └── prediction.py
  ├── app.py              # Streamlit frontend
  ├── run.py              # Script to launch both frontend and backend
  ├── requirements.txt    # Python dependencies
  └── README.md           # Project documentation
```

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for improvements, bug fixes, or new features.

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgments

- Built with Streamlit and FastAPI
- Powered by AI/ML models (scikit-learn)
- Integrated with UniProt, AlphaFold, and public bioinformatics databases 
=======
# MUTALYZE
>>>>>>> 411608e374ed4c6ee7250f1e3b28a93c4e11fcd7
