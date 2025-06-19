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
from datetime import datetime

class MutationPredictor:
    def __init__(self):
        self.impact_levels = {
            'HIGH': ['Pathogenic', 'Likely Pathogenic'],
            'MODERATE': ['Uncertain Significance'],
            'LOW': ['Likely Benign', 'Benign']
        }
        self.uniprot_api = "https://rest.uniprot.org/uniprotkb/search"
        self.alphafold_api = "https://alphafold.ebi.ac.uk/api/prediction/"
        self.interpro_api = "https://rest.uniprot.org/uniprotkb/search"
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
            # Query UniProt API with a simpler query first
            params = {
                'query': f'gene:{gene_name} AND reviewed:true',
                'format': 'json',
                'fields': 'accession,id,sequence,protein_name,gene_names,length,ft_domain,ft_motif,ft_region,ft_repeat,ft_zn_fing,ft_dna_bind,ft_act_site,ft_binding,ft_site,ft_disulfid,ft_crosslnk,ft_mutagen,ft_variant,cc_function,cc_subcellular_location'
            }
            
            # Add logging for the exact API URL and parameters
            print(f"UniProt API URL: {self.uniprot_api}")
            print(f"UniProt API Params: {params}")
            
            response = requests.get(self.uniprot_api, params=params)
            response.raise_for_status()
            
            data = response.json()
            
            # Log the full response data structure
            print(f"UniProt API Response Data Type: {type(data)}")
            print(f"UniProt API Response Data: {json.dumps(data, indent=2)}")

            if not data.get('results'):
                print(f"UniProt API returned no results for gene: {gene_name}")
                return None
                
            result = data['results'][0]
            
            # Log the result object type and keys
            print(f"UniProt API Result Object Type: {type(result)}")
            print(f"UniProt API Result Object Keys: {result.keys()}")

            # Extract protein features
            features = {}
            if 'features' in result and isinstance(result['features'], list):
                for feature in result['features']:
                    # Ensure the feature is a dictionary before processing
                    if isinstance(feature, dict):
                         # Safely get the type and add the feature if type exists
                         feature_type = feature.get('type')
                         if feature_type:
                              features[feature_type] = features.get(feature_type, []) + [feature]
            
            # Filter for relevant feature types
            relevant_features = {k: v for k, v in features.items() if k in ['Domain', 'Motif', 'Region', 'Repeat', 'Zn-finger', 'DNA binding', 'Active site', 'Binding site', 'Site', 'Disulfide bond', 'Cross-link', 'Mutagenesis', 'Variant']}

            # Safely extract nested fields
            sequence_value = result.get('sequence', {}).get('value', '')
            print(f"Extracted sequence_value: {sequence_value[:50]}...") # Log extracted value

            uniprot_id = result.get('primaryAccession', '')
            print(f"Extracted uniprot_id: {uniprot_id}") # Log extracted value
            
            protein_description = result.get('proteinDescription', {})
            print(f"Extracted protein_description type: {type(protein_description)}") # Log type

            recommended_name = protein_description.get('recommendedName', {})
            print(f"Extracted recommended_name type: {type(recommended_name)}") # Log type

            full_name = recommended_name.get('fullName', {})
            print(f"Extracted full_name type: {type(full_name)}") # Log type

            protein_name = full_name.get('value', '')
            print(f"Extracted protein_name (from recommended): {protein_name}") # Log value
            
            # Fallback to submission names if recommended name is not available
            if not protein_name:
                 submission_names = protein_description.get('submissionNames', [])
                 print(f"Extracted submission_names type: {type(submission_names)}") # Log type
                 if submission_names and isinstance(submission_names, list) and isinstance(submission_names[0], dict):
                      protein_name = submission_names[0].get('fullName', '')
                      print(f"Extracted protein_name (from submission): {protein_name}") # Log value

            genes = result.get('genes', [])
            print(f"Extracted genes type: {type(genes)}") # Log type
            gene_names = ''
            if genes and isinstance(genes, list) and isinstance(genes[0], dict):
                 gene_name_obj = genes[0].get('geneName', {})
                 print(f"Extracted gene_name_obj type: {type(gene_name_obj)}") # Log type
                 gene_names = gene_name_obj.get('value', '')
                 print(f"Extracted gene_names: {gene_names}") # Log value

            sequence_length = result.get('sequence', {}).get('length', 0)
            print(f"Extracted sequence_length: {sequence_length}") # Log value
            
            # Safely extract function and subcellular location from comments
            comments = result.get('comments', []) # Get comments, default to list
            print(f"Extracted comments type: {type(comments)}") # Log type

            function_comments = []
            subcellular_location_comments = []

            # Iterate through comments and extract FUNCTION and SUBCELLULAR LOCATION
            if isinstance(comments, list):
                 for comment in comments:
                      # Log comment type and content
                      print(f"Processing comment type: {type(comment)}, content keys: {comment.keys() if isinstance(comment, dict) else 'N/A'}")

                      if isinstance(comment, dict) and comment.get('commentType') == 'FUNCTION':
                           if 'texts' in comment and isinstance(comment['texts'], list):
                                function_comments = [text.get('value', '') for text in comment['texts'] if isinstance(text, dict)]
                      elif isinstance(comment, dict) and comment.get('commentType') == 'SUBCELLULAR LOCATION':
                           if 'locations' in comment and isinstance(comment['locations'], list):
                                # Extract location string from description within location
                                subcellular_location_comments = [location.get('location', {}).get('value', '') for location in comment['locations'] if isinstance(location, dict)]

            print(f"Final function_comments: {function_comments}") # Log final extracted comments
            print(f"Final subcellular_location_comments: {subcellular_location_comments}") # Log final extracted comments

            return {
                'sequence': sequence_value,
                'uniprot_id': uniprot_id,
                'protein_name': protein_name,
                'gene_names': gene_names,
                'length': sequence_length,
                'features': relevant_features, # Store only relevant features
                'function': function_comments, # Use the safely extracted comments
                'subcellular_location': subcellular_location_comments # Use the safely extracted comments
            }
        except requests.exceptions.RequestException as e:
            print(f"Error fetching protein sequence from UniProt API: {e}")
            if e.response is not None:
                print(f"UniProt API response status code: {e.response.status_code}")
                print(f"UniProt API response body: {e.response.text}")
            return None
        except Exception as e:
            print(f"An unexpected error occurred while fetching protein sequence: {e}")
            return None

    def fetch_protein_structure(self, uniprot_id: str) -> Optional[str]:
        """Fetch protein structure from AlphaFold"""
        try:
            response = requests.get(f"{self.alphafold_api}{uniprot_id}")
            response.raise_for_status()
            # Assuming AlphaFold API returns a list and the first element contains the pdbUrl
            data = response.json()
            if isinstance(data, list) and data:
                 return data[0].get('pdbUrl')
            return None
        except Exception as e:
            print(f"Error fetching protein structure: {str(e)}")
            return None

    def calculate_conservation_score(self, mutation: str, sequence: str) -> float:
        """Calculate conservation score using multiple sequence alignment"""
        try:
            # Query NCBI BLAST for similar sequences
            position = int(mutation[1:-1])
            window = 10  # Look at surrounding amino acids
            
            # Calculate local conservation
            start = max(0, position - window)
            end = min(len(sequence), position + window)
            local_region = sequence[start:end]
            
            # Query UniProt for similar sequences
            params = {
                'query': f'sequence:"{local_region}" AND reviewed:true',
                'format': 'json',
                'fields': 'accession,sequence',
                'size': 100
            }
            response = requests.get(self.uniprot_api, params=params)
            response.raise_for_status()
            
            data = response.json()
            if not data['results']:
                print(f"UniProt API returned no results for sequence: {local_region}")
                return 0.5
            
            # Calculate conservation score based on sequence alignment (simplified)
            # This is a basic approach; real MSA tools are more complex
            aligned_sequences = [result['sequence']['value'] for result in data['results']]
            conservation_scores = []
            
            if not aligned_sequences:
                return 0.5 # Return default if no similar sequences found

            try:
                 for i in range(len(local_region)):
                    column = [seq[position - len(local_region) + i] if position - len(local_region) + i < len(seq) and position - len(local_region) + i >= 0 else '-' for seq in aligned_sequences]
                    if column:
                        most_common = max(set(column), key=column.count)
                        conservation = column.count(most_common) / len(column)
                        conservation_scores.append(conservation)

                 if not conservation_scores:
                    return 0.5 # Return default if no scores calculated

                 return sum(conservation_scores) / len(conservation_scores)

            except Exception as e:
                 print(f"Error calculating conservation score during alignment processing: {str(e)}")
                 return 0.5
            
        except requests.exceptions.RequestException as e:
            print(f"Error fetching similar sequences from UniProt API: {e}")
            if e.response is not None:
                print(f"UniProt API response status code: {e.response.status_code}")
                print(f"UniProt API response body: {e.response.text}")
            return 0.5
        except Exception as e:
            print(f"An unexpected error occurred while calculating conservation score: {e}")
            return 0.5

    def analyze_mutation(self, gene_name: str, mutation: str) -> Dict:
        """Analyze mutation and predict impact using rule-based logic"""
        # Extract mutation details
        original_aa = mutation[0]
        position = int(mutation[1:-1])
        new_aa = mutation[-1]
        
        # Get protein sequence and metadata
        protein_data = self.fetch_protein_sequence(gene_name)
        if not protein_data:
            # Return a complete response with default values when protein data can't be fetched
            return {
                'gene': gene_name,
                'mutation': mutation,
                'impact': 'Unknown',
                'confidence': 0.0,
                'conservation_score': 0.0,
                'protein_domain': 'Unknown',
                'clinical_significance': 'Could not fetch protein data',
                'protein_name': f"{gene_name} protein",
                'uniprot_id': 'Unknown',
                'structure_url': None,
                'mutation_details': {
                    'position': position,
                    'original': original_aa,
                    'mutated': new_aa,
                    'type': 'Missense'
                },
                'protein_features': {}, # Provide empty dict for features
                'function': [],
                'subcellular_location': [],
                'error': 'Could not fetch protein sequence'
            }
        
        # Calculate various scores
        conservation_score = self.calculate_conservation_score(mutation, protein_data['sequence'])
        
        # Get protein structure if available
        structure_url = self.fetch_protein_structure(protein_data['uniprot_id'])
        
        # Determine protein domain
        domain = self._get_protein_domain(protein_data['features'], position)
        
        # Get clinical significance based on features and conservation
        clinical_significance = self._get_clinical_significance(protein_data['features'], position, conservation_score, original_aa, new_aa)
        
        # Calculate rule-based confidence score (simplified)
        # This is a basic example; a real confidence score would be more complex
        confidence = conservation_score # Use conservation score as a proxy for confidence

        # Safely extract protein name and gene names
        protein_name = protein_data.get('protein_name', f"{gene_name} protein")
        gene_names = protein_data.get('gene_names', '')

        # Safely extract function and subcellular_location, handling potential list format
        function_comments = []
        if isinstance(protein_data.get('function'), list):
             function_comments = protein_data.get('function')
        elif isinstance(protein_data.get('function'), dict) and 'texts' in protein_data.get('function'):
             function_comments = [text['value'] for text in protein_data['function']['texts']]
        
        subcellular_location_comments = []
        if isinstance(protein_data.get('subcellular_location'), list):
             subcellular_location_comments = protein_data.get('subcellular_location')
        elif isinstance(protein_data.get('subcellular_location'), dict) and 'texts' in protein_data.get('subcellular_location'):
             subcellular_location_comments = [text['value'] for text in protein_data['subcellular_location']['texts']]

        return {
            'gene': gene_name,
            'mutation': mutation,
            'impact': clinical_significance,
            'confidence': float(confidence), # Ensure confidence is a float
            'conservation_score': conservation_score,
            'protein_domain': domain,
            'clinical_significance': clinical_significance,
            'protein_name': protein_name, # Use the safely extracted protein name
            'uniprot_id': protein_data['uniprot_id'],
            'structure_url': structure_url,
            'mutation_details': {
                'position': position,
                'original': original_aa,
                'mutated': new_aa,
                'type': 'Missense'
            },
            'protein_features': protein_data['features'], # Include features in the response
            'function': function_comments, # Use the safely extracted comments
            'subcellular_location': subcellular_location_comments # Use the safely extracted comments
        }
    
    def _get_protein_domain(self, features: Dict, position: int) -> str:
        """Get protein domain information from features"""
        try:
            for feature_type, feature_list in features.items():
                for feature in feature_list:
                    # Check if the feature has a location and description
                    if 'location' in feature and 'description' in feature:
                        try:
                            start = int(feature['location']['start']['value'])
                            end = int(feature['location']['end']['value'])
                            # If the mutation position is within the feature's location, return its description
                            if start <= position <= end:
                                return f"{feature_type}: {feature['description']}"
                        except (ValueError, TypeError) as e:
                             print(f"Error parsing feature location for {feature_type}: {e}")
                             continue # Skip to the next feature if location parsing fails
            # If no specific feature is found at the position, return a general message
            return "No specific domain or feature found at this position"
        except Exception as e:
            print(f"Error getting protein domain: {str(e)}")
            return "Error retrieving protein domain information"
    
    def _get_clinical_significance(self, features: Dict, position: int, conservation_score: float, original_aa: str, new_aa: str) -> str:
        """Get clinical significance based on features and conservation score (rule-based)"""
        try:
            # Rule 1: Check for known pathogenic or disease-associated variants in features
            for feature_type, feature_list in features.items():
                for feature in feature_list:
                    if 'location' in feature and 'description' in feature:
                        try:
                            start = int(feature['location']['start']['value'])
                            end = int(feature['location']['end']['value'])
                            if start <= position <= end:
                                description = feature['description'].lower()
                                if 'pathogenic' in description or 'disease' in description or 'clinical significance' in description:
                                    # Prioritize direct clinical annotations
                                    return f"Likely Pathogenic (annotated {feature_type}: {feature['description']})"
                        except (ValueError, TypeError) as e:
                             print(f"Error parsing feature location for significance: {e}")
                             continue

            # Rule 2: Consider mutations in highly conserved functional sites/domains
            if conservation_score > 0.8:
                for feature_type in ['Active site', 'Binding site', 'Zn-finger', 'DNA binding', 'Domain']:
                     if feature_type in features:
                           for feature in features[feature_type]:
                                if 'location' in feature:
                                    try:
                                         start = int(feature['location']['start']['value'])
                                         end = int(feature['location']['end']['value'])
                                         if start <= position <= end:
                                              return f"Likely Pathogenic (high conservation in {feature_type})"
                                    except (ValueError, TypeError) as e:
                                         print(f"Error parsing feature location for conservation rule: {e}")
                                         continue

            # Rule 3: Consider the type of amino acid change and conservation
            amino_acid_properties = {
                'A': 'Nonpolar', 'V': 'Nonpolar', 'I': 'Nonpolar', 'L': 'Nonpolar', 'M': 'Nonpolar', 'F': 'Nonpolar', 'W': 'Nonpolar', 'P': 'Nonpolar',
                'G': 'Nonpolar',
                'S': 'Polar', 'T': 'Polar', 'C': 'Polar', 'Y': 'Polar', 'N': 'Polar', 'Q': 'Polar',
                'D': 'Acidic', 'E': 'Acidic',
                'K': 'Basic', 'R': 'Basic', 'H': 'Basic'
            }

            original_prop = amino_acid_properties.get(original_aa, 'Unknown')
            new_prop = amino_acid_properties.get(new_aa, 'Unknown')

            # Large change in amino acid property in a moderately conserved region
            if original_prop != new_prop and conservation_score > 0.6:
                 return f"Uncertain Significance (change from {original_prop} to {new_prop} in conserved region)"

            # Rule 4: Default based on conservation score if no specific features or significant property changes
            if conservation_score > 0.7:
                return "Uncertain Significance (moderately conserved region)"
            elif conservation_score > 0.4:
                return "Likely Benign (low to moderate conservation)"
            else:
                return "Benign (low conservation)"

        except Exception as e:
            print(f"Error getting clinical significance: {str(e)}")
            return "Error determining clinical significance" 