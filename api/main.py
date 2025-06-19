from fastapi import FastAPI, HTTPException, Depends
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional, Dict
import uvicorn
from datetime import datetime
import os
from dotenv import load_dotenv
from utils.prediction import MutationPredictor

# Load environment variables
load_dotenv()

# Initialize FastAPI app
app = FastAPI(
    title="Mutalyze API",
    description="API for gene mutation impact prediction",
    version="1.0.0"
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize prediction model
predictor = MutationPredictor()

# In-memory storage
predictions_db: List[Dict] = []

# Models
class MutationRequest(BaseModel):
    gene_name: str
    mutation: str

class BatchMutationRequest(BaseModel):
    mutations: List[MutationRequest]

class MutationResponse(BaseModel):
    impact: str
    confidence: float
    conservation_score: float
    protein_domain: str
    clinical_significance: str
    protein_name: Optional[str]
    uniprot_id: Optional[str]
    structure_url: Optional[str]
    mutation_details: dict
    timestamp: datetime

# Routes
@app.get("/")
async def root():
    return {"message": "Welcome to Mutalyze API"}

@app.post("/predict", response_model=MutationResponse)
async def predict_mutation(request: MutationRequest):
    try:
        # Get prediction
        result = predictor.analyze_mutation(request.gene_name, request.mutation)
        
        # Add timestamp
        result['timestamp'] = datetime.now()
        
        # Store in memory
        predictions_db.append(result)
        
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/batch-predict", response_model=List[MutationResponse])
async def batch_predict(request: BatchMutationRequest):
    try:
        results = []
        for mutation in request.mutations:
            result = predictor.analyze_mutation(mutation.gene_name, mutation.mutation)
            result['timestamp'] = datetime.now()
            results.append(result)
            
            # Store in memory
            predictions_db.append(result)
            
        return results
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/history")
async def get_history(limit: int = 10):
    try:
        # Return last N predictions
        return sorted(predictions_db, key=lambda x: x['timestamp'], reverse=True)[:limit]
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/protein/{uniprot_id}")
async def get_protein_info(uniprot_id: str):
    try:
        # Get protein sequence and metadata
        protein_data = predictor.fetch_protein_sequence(uniprot_id)
        if not protein_data:
            raise HTTPException(status_code=404, detail="Protein not found")
        return protein_data
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True) 