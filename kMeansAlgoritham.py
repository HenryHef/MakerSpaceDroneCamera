import random
import numpy as np

def init_board_gauss(N, k):
    n = float(N)/k
    X = []
    for i in range(k):
        c = (random.uniform(-1, 1), random.uniform(-1, 1))
        s = random.uniform(0.05,0.5)
        x = []
        while len(x) < n:
            a, b = np.array([np.random.normal(c[0], s), np.random.normal(c[1], s)])
            # Continue drawing points from the distribution in the range [-1,1]
            if abs(a) < 1 and abs(b) < 1:
                x.append([a,b])
        X.extend(x)
    X = np.array(X)[:N]
    return X

def cluster_points(X, mu):
    print "clustering "+str(len(X))
    clusters  = {}
    for x in X:
        bestmukey = min([(i[0], np.linalg.norm([x[j]-mu[i[0]][j] for j in range(len(x))])) for i in enumerate(mu)], key=lambda t:t[1])[0]
        try:
            clusters[bestmukey].append(x)
        except KeyError:
            clusters[bestmukey] = [x]
    print "end clustering"
    return clusters
 
def reevaluate_centers(mu, clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        newmu.append(np.mean(clusters[k], axis = 0))
    return newmu

def has_converged(mu, oldmu):
	return set([tuple(a) for a in mu]) == set([tuple(a) for a in oldmu])


def find_centers(X, K):
    N=20 # helps with points hovering between two clusters
	# Initialize Data=X to K random centers
    oldmu = random.sample(X, K*N)
    mu = random.sample(X, K*N)
    count = 0
    while not has_converged(mu, oldmu) or count<=3:
        oldmu = mu
        # Assign all points in X to clusters
        clusters = cluster_points(X, mu)
        # Reevaluate centers
        mu = reevaluate_centers(oldmu, clusters)
        if count==3 or has_converged(mu, oldmu):#hack to avoid only centers in the same zone
            muC = mu
            mu=[]
            for elem in muC:
                if all((e[0]-elem[0])**2+(e[1]-elem[1])**2>60**2 for e in mu):
                    mu.append(elem)
            count=3
            print mu
            mu=random.sample(mu,K)
        count=count+1
    return(mu, clusters)
