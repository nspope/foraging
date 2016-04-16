modelcode <- list()
## model 1: without floral assessment
	modelcode[[1]] <- '
		data {
			int<lower=1> C;
			int<lower=1> K;
			matrix[K, 2] trap;
			matrix[C, K] y;
			real lowerbound;
			real upperbound;
			vector[K] floral;
			real<lower=0> priorVa;
			real<lower=0> priorCo;
		}
		parameters {
			real<lower=0> beta;
			real<lower=0> sigma;
			real<lower=0> tau;
			real mu;
			vector[K] eps;
			vector[C] zeta;
			matrix<lower=lowerbound, upper=upperbound>[C,2] delta;
		}
		transformed parameters {
			real<lower=0> tau_sqrt;
			real<lower=0> sigma_sqrt;
			vector[C] zeta_scale;
			vector[K] eps_scale;

			tau_sqrt <- sqrt(tau);
			sigma_sqrt <- sqrt(sigma);
			eps_scale <- eps*sigma_sqrt;
			zeta_scale <- zeta*tau_sqrt;
		}
		model {
			// temporary declarations
			matrix[C,K] dis;
			matrix[C,K] lambda;

			// priors
			sigma ~ normal(0, priorVa);			
			tau ~ normal(0, priorVa);
			beta ~ normal(0, priorCo);
			mu ~ normal(0, priorCo);

			// random effects for traps
			eps ~ normal(0, 1);
			zeta ~ normal(0, 1);

			// calculate intensity
			for(k in 1:K){	
				for(i in 1:C){
					dis[i, k] <- distance(delta[i], trap[k]);
					lambda[i, k] <- -beta*dis[i, k] + mu + zeta_scale[i] + eps_scale[k];
				}
			}

			// log likelihood
			increment_log_prob( -exp(lambda) );
			increment_log_prob( y .* lambda );
		}
'


## model 2: local floral assessment	
	modelcode[[2]] <- '
		data {
			int<lower=1> C;
			int<lower=1> K;
			matrix[K, 2] trap;
			matrix[C, K] y;
			real lowerbound;
			real upperbound;
			vector[K] floral;
			real<lower=0> priorVa;
			real<lower=0> priorCo;
		}
		parameters {
			real<lower=0> beta;
			real<lower=0> sigma;
			real<lower=0> tau;
			real theta;
			real mu;
			vector[K] eps;
			vector[C] zeta;
			matrix<lower=lowerbound, upper=upperbound>[C,2] delta;
		}
		transformed parameters {
			real<lower=0> tau_sqrt;
			real<lower=0> sigma_sqrt;
			vector[C] zeta_scale;
			vector[K] eps_scale;

			tau_sqrt <- sqrt(tau);
			sigma_sqrt <- sqrt(sigma);
			eps_scale <- eps*sigma_sqrt;
			zeta_scale <- zeta*tau_sqrt;
		}
		model {
			// temporary declarations
			matrix[C,K] dis;
			matrix[C,K] lambda;

			// priors
			sigma ~ normal(0, priorVa);			
			tau ~ normal(0, priorVa);
			beta ~ normal(0, priorCo);
			mu ~ normal(0, priorCo);
			theta ~ normal(0, priorCo);

			// random effects for traps
			eps ~ normal(0, 1);
			zeta ~ normal(0, 1);

			// calculate intensity
			for(k in 1:K){	
				for(i in 1:C){
					dis[i, k] <- distance(delta[i], trap[k]);
					lambda[i, k] <- -beta*dis[i, k] + theta*floral[k] + mu + zeta_scale[i] + eps_scale[k];
				}
			}

			// log likelihood
			increment_log_prob( -exp(lambda) );
			increment_log_prob( y .* lambda );
		}
'

