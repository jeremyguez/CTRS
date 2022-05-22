
function (void)outputValues(numeric alpha, numeric K_used, numeric vecM, numeric vecE, numeric vecMG, numeric vecEG, numeric vecMF, numeric vecEF, numeric vecHistEnf)
{
	
	cat(sim.generation + " " + alpha + " " + K_used + " " + 
	var(vecM) + " " + var(vecE) + " " + cov(vecM,vecE) + " " +
	cor(vecMF,vecEF) + " " + cor(vecMG,vecEG) + " " + cor(vecM,vecE) + " ");
	
	for (i in seqLen(size(vecHistEnf)))
	{
		cat(asString(vecHistEnf[i]) + " ");
	}
	catn();
	
}

// Each individual has a property "tagF" which contains its sibshipp size
// apart from first generation where those tags are randomly drawn

initialize() {
	
	initializeSLiMModelType("nonWF"); // non WF simulation
	initializeSLiMOptions(keepPedigrees=T); // keep pedigree to prevent inbreeding		
	initializeSex("A"); // autosome
	defineConstant("K", sizeFirst_VAL);	// K = Pop size most ancient
	defineConstant("K2", sizeChange_VAL);	// K2 = Pop size after one change (if there is a change)
	defineConstant("K3", size2Change_VAL);	// K3 = Pop size after a seconde demographic change (if there is a change)
	defineConstant("GENCHANGE", genChange_VAL); // Generation when first demography change occurs
	defineConstant("GENCHANGE2", genChange2_VAL); // Generation when second demography change occurs
	defineConstant("ENR", c(steps_VAL)); // steps of outputs
	defineConstant("ALPHA", alpha_VAL); // alpha
	defineConstant("DISTRIB", "distrib_Bash"); // distribution type (poisson or geometric)
	defineConstant("BURNING", burnin_VAL); // length of burning
	defineConstant("EVENT", "event_Bash"); // type of demographic event
	defineConstant("FACTOR", factor_VAL); // factor for exponential event
	defineConstant("EVENT2", "event2_Bash"); // type of demographic event 2
	defineConstant("FACTOR2", factor2_VAL); // factor for exponential event 2
	


	// neutral mutations
	initializeMutationType("m1", 0.5, "f", 0.0);
	m1.convertToSubstitution = T;
	initializeTreeSeq(); // keeping trees
	
	initializeGenomicElementType("g1", m1, 1.0);
	initializeGenomicElement(g1, 0, L_VAL); // genome size
	initializeMutationRate(0); // No mutations, added on the trees later on if needed
	initializeRecombinationRate(rho_VAL); // recombination rate
}

1 early() {
	
	if(GENCHANGE < 1 & GENCHANGE2 > 1)
	{
		K_used = K2;
	}
	else if(GENCHANGE2 < 1)
	{
		K_used = K3;
	}
	
	sim.addSubpop("p1", K_used);
	
	p1.individuals.tagF = rep(1, p1.individualCount);
	
	p1.individuals.x = rep(0, p1.individualCount);
	
	
}


reproduction(NULL, NULL) 
{
	
	self.active = 0;  // one call of "reproduction" per generation
	
	inds = p1.individuals;
	l = length(inds);
	
	
	if(EVENT == "Sudden")
	{
		if (sim.generation >= GENCHANGE & sim.generation < GENCHANGE2) // if we are after the change of pop size
		{
			K_used = K2;
		}
		else if (sim.generation >= GENCHANGE2)
		{
			K_used = K3;
		}
	}
	else if (EVENT == "Exponential")
	{
		if (sim.generation >= GENCHANGE & sim.generation < GENCHANGE2) // if we are after the change of pop size
		{
			if(K2 > K & l < K2) // if it's an expansion AND K_used still did not reached K2
			{
				
				K_used = l * FACTOR;
				if(K_used >= K2)
				{
					K_used = K2;
				}
			
			}
			else if(K2 < K & l > K2) // if it's a contraction AND K_used still did not reached K2
			{
				K_used = l * FACTOR;
				if(K_used >= K2)
				{
					K_used = K2;
				}
				
			}
			else
			{
				K_used = l;
			}	
			
		}
		else if (sim.generation >= GENCHANGE2)
		{
			if(K3 > K2 & l < K3) // if it's an expansion AND K_used still did not reached  K2
			{
				
				K_used = l * FACTOR2;
				if(K_used >= K3)
				{
					K_used = K3;
				}
			
			}
			else if(K3 < K2 & l > K3) // if it's a contraction AND K_used still did not reached K2
			{
				K_used = l * FACTOR2;
				if(K_used >= K3)
				{
					K_used = K3;
				}
				
			}
			else
			{
				K_used = l;
			}
		}
	}
	

	if(isFloat(K_used))
	{
		K_used = asInteger(ceil(K_used));
	}



	//######## Section 1 : creating couples ########//
	
	
	
    indsF = inds[inds.sex=='F'];
    indsM = inds[inds.sex=='M'];
    
    nb_F = length(indsF); // nb of women
	nb_M = length(indsM); // nb of men
	
	
    // alpha computation
    
    debut = deb_Bash; 
	fin = fin_Bash;
    
    if (sim.generation > debut & sim.generation <= fin)
	{
		alpha = ALPHA;
	}
	else
	{
		alpha = 0;
	}	


	if (nb_F < nb_M) 
	{
		nb_couples = nb_F;
	}
	else
	{
		nb_couples = nb_M;
	}
	
	
	indsF = p1.sampleIndividuals(nb_couples, sex='F');
	indsM = p1.sampleIndividuals(nb_couples, sex='M');
	
	
	
	//### Section  2 : computing probability of reproduction for each couple ###//
	
	if (DISTRIB == "P") // if distribution is Poisson like
	{
		distrib = rep(1,nb_couples);
	}
	else
	{
		distrib = rgamma(nb_couples, 1, 1);
	}
	
	proba = distrib * (((indsM.tagF + indsF.tagF)/2)^(alpha));
	
	
	//### Section 3 : creating K children ###//
	
	
	
	listeEnfants=c();
	sexe = rep(c("M","F"), asInteger(ceil(K_used/2)));
	
	
	parent_indices = sample(x = 0:(length(indsM)-1), size = K_used, replace = T, weights = proba);
	
	listeEnfants = sapply(0:(K_used-1), "subpop.addCrossed(indsF[parent_indices[applyValue]], indsM[parent_indices[applyValue]], sex=sexe[applyValue]);");
	
	indsM.tag = indsM.pedigreeID; // keeps pedigree in the tag
	
	for(i in seqLen(K_used))
	{
		parent_index = parent_indices[i];
		
		indsM[parent_index].x = indsM[parent_index].x + 1.0; // incrementing x (nb of children) for the father
	}
	
	indsF.x = indsM.x;
	
	
	listeEnfants.x = 0.0; // chlidren start with 0
	
	
	
		
	listeEnfants.tagF = indsM[parent_indices].x;

	if(sim.generation >= 2 & any(sim.generation == ENR)) 
	{
		
		vecM = c();
		vecF = c();
		vecE = c();
		for (i in seqLen(size(inds)))
		{
			vecM = c(vecM, inds[i].tagF); 
			vecE = c(vecE, inds[i].x);		
		}
		
		vecMG = c();
		vecFG = c();
		vecEG = c();
		
		for (i in seqLen(size(indsM)))
		{
			vecMG = c(vecMG, indsM[i].tagF);
			vecEG = c(vecEG, indsM[i].x);	
		}
		
		
		
		vecMF = c();
		vecFF = c();
		vecEF = c();
		
		
		for (i in seqLen(size(indsF)))
		{
			vecMF = c(vecMF, indsF[i].tagF);
			vecEF = c(vecEF, indsF[i].x);			
		}
		
		
		nb = nb_Bash;
		vecHistEnf = rep(0,nb+1);
		
		for (i in seqLen(nb_couples))
		{
			
			val = asInteger(indsM[i].x[0]);
			if (val < nb)
			{
				vecHistEnf[val] = vecHistEnf[val] + 1;
			}
			else
			{
				vecHistEnf[nb] = vecHistEnf[nb] + 1;
			}
		}
			
			
		vecHistEnf2 = vecHistEnf / nb_couples;
		
		// Outputting values for analysis
		outputValues(alpha, K_used, vecM, vecE, vecMG, vecEG, vecMF, vecEF, vecHistEnf2); 
	}
    
}



early() {
	
	inds = p1.individuals;
	adults = inds[inds.age > 0];
	
	// non-overlapping generations : killing adults
	adults.fitnessScaling = 0.0;
	
	
}


