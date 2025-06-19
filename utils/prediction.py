import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import requests
from typing import Dict, Tuple, Optional
import json
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
import joblib
import os

class MutationPredictor:
    def __init__(self):
        self.impact_levels = {
            'HIGH': ['Pathogenic', 'Likely Pathogenic'],
            'MODERATE': ['Uncertain Significance'],
            'LOW': ['Likely Benign', 'Benign']
        }
        self.uniprot_api = "https://rest.uniprot.org/uniprotkb/search"
        self.alphafold_api = "https://alphafold.ebi.ac.uk/api/prediction/"
        self.model_path = "models/mutation_predictor.joblib"
        self._load_or_train_model()
        
    def _load_or_train_model(self):
        """Load or train the mutation prediction model"""
        if os.path.exists(self.model_path):
            self.model = joblib.load(self.model_path)
        else:
            # Train a simple model (in production, use real training data)
            self.model = RandomForestClassifier(n_estimators=100)
            # Mock training data
            X = np.random.rand(100, 5)  # Features: conservation, position, etc.
            y = np.random.choice(['Pathogenic', 'Benign'], 100)
            self.model.fit(X, y)
            os.makedirs("models", exist_ok=True)
            joblib.dump(self.model, self.model_path)
        
    def fetch_protein_sequence(self, gene_name: str) -> Optional[Dict]:
        """Fetch protein sequence and metadata from UniProt API"""
        try:
            # Query UniProt API
            params = {
                'query': f'gene:{gene_name} AND reviewed:yes',
                'format': 'json',
                'fields': 'sequence,id,protein_name,gene_names,length'
            }
            response = requests.get(self.uniprot_api, params=params)
            response.raise_for_status()
            
            data = response.json()
            if not data['results']:
                return None
                
            result = data['results'][0]
            return {
                'sequence': result['sequence']['value'],
                'uniprot_id': result['primaryAccession'],
                'protein_name': result['proteinDescription']['recommendedName']['fullName']['value'],
                'gene_names': result['genes'][0]['geneName']['value'],
                'length': result['sequence']['length']
            }
        except Exception as e:
            print(f"Error fetching protein sequence: {str(e)}")
            return None

    def fetch_protein_structure(self, uniprot_id: str) -> Optional[str]:
        """Fetch protein structure from AlphaFold"""
        try:
            response = requests.get(f"{self.alphafold_api}{uniprot_id}")
            response.raise_for_status()
            return response.json()['pdbUrl']
        except Exception as e:
            print(f"Error fetching protein structure: {str(e)}")
            return None

    def calculate_conservation_score(self, mutation: str, sequence: str) -> float:
        """Calculate conservation score using multiple sequence alignment"""
        try:
            # Query NCBI BLAST for similar sequences
            # This is a simplified version - in production, use actual BLAST API
            position = int(mutation[1:-1])
            window = 10  # Look at surrounding amino acids
            
            # Calculate local conservation
            start = max(0, position - window)
            end = min(len(sequence), position + window)
            local_region = sequence[start:end]
            
            # Simple conservation score based on amino acid properties
            amino_acid_properties = {
                'A': 0.8, 'V': 0.8, 'I': 0.8, 'L': 0.8,  # Hydrophobic
                'F': 0.9, 'W': 0.9, 'Y': 0.9,  # Aromatic
                'S': 0.7, 'T': 0.7, 'N': 0.7, 'Q': 0.7,  # Polar
                'D': 0.6, 'E': 0.6,  # Acidic
                'K': 0.6, 'R': 0.6, 'H': 0.6,  # Basic
                'C': 0.9, 'G': 0.5, 'P': 0.5, 'M': 0.8  # Special cases
            }
            
            # Calculate weighted conservation score
            scores = [amino_acid_properties.get(aa, 0.5) for aa in local_region]
            return sum(scores) / len(scores)
            
        except Exception as e:
            print(f"Error calculating conservation score: {str(e)}")
            return 0.5

    def analyze_mutation(self, gene_name: str, mutation: str) -> Dict:
        """Analyze mutation and predict impact using ML model"""
        # Extract mutation details
        original_aa = mutation[0]
        position = int(mutation[1:-1])
        new_aa = mutation[-1]
        
        # Get protein sequence and metadata
        protein_data = self.fetch_protein_sequence(gene_name)
        if not protein_data:
            return {
                'impact': 'Unknown',
                'confidence': 0.0,
                'error': 'Could not fetch protein sequence'
            }
        
        # Calculate various scores
        conservation_score = self.calculate_conservation_score(mutation, protein_data['sequence'])
        
        # Prepare features for ML model
        features = np.array([[
            conservation_score,
            position / protein_data['length'],  # Normalized position
            len(protein_data['sequence']),  # Protein length
            amino_acid_properties.get(original_aa, 0.5),  # Original AA property
            amino_acid_properties.get(new_aa, 0.5)  # New AA property
        ]])
        
        # Get prediction from model
        impact_probs = self.model.predict_proba(features)[0]
        impact_idx = np.argmax(impact_probs)
        confidence = impact_probs[impact_idx]
        
        # Map prediction to impact level
        impact = self.model.classes_[impact_idx]
        
        # Get protein structure if available
        structure_url = self.fetch_protein_structure(protein_data['uniprot_id'])
        
        return {
            'impact': impact,
            'confidence': float(confidence),
            'conservation_score': conservation_score,
            'protein_domain': self._get_protein_domain(gene_name, position),
            'clinical_significance': self._get_clinical_significance(impact),
            'protein_name': protein_data['protein_name'],
            'uniprot_id': protein_data['uniprot_id'],
            'structure_url': structure_url,
            'mutation_details': {
                'original_aa': original_aa,
                'position': position,
                'new_aa': new_aa
            }
        }
    
    def _get_protein_domain(self, gene_name: str, position: int) -> str:
        """Get protein domain information from InterPro"""
        try:
            # Query InterPro API for domain information
            # This is a mock implementation - in production, use actual InterPro API
            domains = {
                'BRCA1': {
                    'RING': (1, 100),
                    'BRCT': (1600, 1863)
                },
                'TP53': {
                    'DNA_binding': (100, 300),
                    'Tetramerization': (325, 355)
                }
            }
            
            if gene_name in domains:
                for domain, (start, end) in domains[gene_name].items():
                    if start <= position <= end:
                        return domain
            return "Unknown domain"
        except Exception as e:
            print(f"Error fetching protein domain: {str(e)}")
            return "Unknown domain"
    
    def _get_clinical_significance(self, impact: str) -> str:
        """Get clinical significance based on impact"""
        if impact in self.impact_levels['HIGH']:
            return "Pathogenic"
        elif impact in self.impact_levels['MODERATE']:
            return "Uncertain"
        else:
            return "Benign" 