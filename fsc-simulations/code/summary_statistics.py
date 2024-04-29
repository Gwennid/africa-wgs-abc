#import data
import numpy as np
import scipy as sp
import warnings
from scipy import stats
#from gwas import missing, complete_cases, is_na
import bisect
from collections import Counter

def spatial_histo_fast(pos_list,count_list,M,dmax=np.inf):
    '''Computes the site frequency spectrum

    Fast version of spatial_histo
    Note: This is the correct implementation of dist
    Author Flora

    Returns :
    - the site frequency spectrum from 1 to M (percentages, sum to 1)
    - the variance of the distance between sites with count i, for all i. This variance needs to be multiplied by the overall proportion of SNPs.
    positions of vector pos are assumed to be SORTED within each chromosome
    '''
    histo=np.zeros(shape=M,dtype='float')
    nb_snp=0
    d = [[] for i in range(1,M+1) ]
    for chro in range(0,len(pos_list)):
        pos= [[] for i in range(1,M+1) ]
        for snp in xrange(pos_list[chro].shape[0]):
            i=count_list[chro][snp]
            try:
                histo[i-1]+=1
                pos[i-1].append(pos_list[chro][snp])
            except IndexError:
                continue

        [d[i-1].append(x) for i in range (1,M+1) for x in np.diff(pos[i-1]) if x<=dmax]

    # for each frequency, compute the std of the distance list, after removing distance longer than dmax
        # if no distances deviation set to 0 (was -1 before)
    dist=np.asarray([np.std(d_at_freq) if len(d_at_freq)>1 else 0.0 for d_at_freq in d])
#    dist=np.asarray([np.std(d_at_freq) if len(d_at_freq)>1 else -1.0 for d_at_freq in d])

    # correct but with np.nan and 0 for len(d)=0 or 1
    # dist=[np.std([x for x in d_at_freq if x<=dmax]) for d_at_freq in d]

    return histo/np.sum(histo),dist


def r2(u,v):
    '''
    returns the r2 value for two haplotype vectors (numpy arrays with alleles coded 1 and 0)
    '''
    fcross=np.mean(u*v)
    fu=np.mean(u)
    fv=np.mean(v)
    return (fcross-fu*fv)**2/(fu*(1-fu)*fv*(1-fv))

def distrib_r2(pos_list,hap_list,interval_list):
    '''
    returns the mean and the variance of r2 for a list of distance intervals.
    pos_list is a list of 1 dim arrays
    hap_list is a list of 2 dim arrays
    interval_list is a list of ordered pairs
    a subset of non overlapping pairs is used for each interval
    '''
    p=len(interval_list)
    moy=-np.ones(shape=p,dtype='float')
    var=-np.ones(shape=p,dtype='float')
    for i in range(0,p):
        r2_list=[]
        dmin=interval_list[i][0]
        dmax=interval_list[i][1]
        # looks for snp pairs with the good distance
        for chro in range(0,len(pos_list)):
            nb_snp=len(pos_list[chro])
            if nb_snp>0:
                i_deb=0
                i_fin=1
                while i_fin<nb_snp:
                    while i_fin < nb_snp and (pos_list[chro][i_fin]-pos_list[chro][i_deb])<dmin:
                        i_fin+=1
                    if i_fin < nb_snp and (pos_list[chro][i_fin]-pos_list[chro][i_deb])<=dmax:
                        # compute r2
                        u_deb=hap_list[chro][:,i_deb]
                        u_fin=hap_list[chro][:,i_fin]
                        r2_list.append(r2(u_deb,u_fin))
                    i_deb=i_fin+1
                    i_fin=i_deb+1
        if len(r2_list) < 2:
            # try a more exhaustive screening of SNP pairs
            r2_list=[]
            dmin=interval_list[i][0]
            dmax=interval_list[i][1]
            for chro in range(0,len(pos_list)):
                nb_snp=len(pos_list[chro])
                if nb_snp>0:
                    i_deb=0
                    i_fin=1
                    while i_fin<nb_snp:
                        while i_fin < nb_snp and (pos_list[chro][i_fin]-pos_list[chro][i_deb])<dmin:
                            i_fin+=1
                        if i_fin < nb_snp and (pos_list[chro][i_fin]-pos_list[chro][i_deb])<=dmax:
                            # compute r2
                            u_deb=hap_list[chro][:,i_deb]
                            u_fin=hap_list[chro][:,i_fin]
                            r2_list.append(r2(u_deb,u_fin))
                        i_deb+=1
                        i_fin=i_deb+1
        # computes the stat
        if len(r2_list) >= 2:
            moy[i]=np.mean(np.array(r2_list,dtype='float'))
            var[i]=np.std(np.array(r2_list,dtype='float'))
    return moy,var

def zyg_r2(u,v):
    '''
    returns the zygotic r2 value for two genotype vectors (numpy arrays with genotypes coded 0, 1 and 2)
    '''
    return (np.corrcoef(u,v)[0,1])**2

def distrib_zyg_r2(pos_list,geno_list,interval_list):
    '''
    returns the mean and the variance of zygotic r2 for a list of distance intervals.
    pos_list is a list of 1 dim arrays
    geno_list is a list of 2 dim arrays
    interval_list is a list of ordered pairs
    a subset of non overlapping pairs is used for each interval
    '''
    warnings.simplefilter("error",RuntimeWarning)
    p=len(interval_list)
    moy=-np.ones(shape=p,dtype='float')
    var=-np.ones(shape=p,dtype='float')
    for i in range(0,p):
        r2_list=[]
        dmin=interval_list[i][0]
        dmax=interval_list[i][1]
        for chro in range(0,len(pos_list)):
            nb_snp=len(pos_list[chro])
            if nb_snp>0:
                i_deb=0
                i_fin=1
                while i_fin<nb_snp:
                    while i_fin < nb_snp and (pos_list[chro][i_fin]-pos_list[chro][i_deb])<dmin:
                        i_fin+=1
                    if i_fin < nb_snp and (pos_list[chro][i_fin]-pos_list[chro][i_deb])<=dmax:
                        # compute r2
                        u_deb=geno_list[chro][:,i_deb]
                        u_fin=geno_list[chro][:,i_fin]
                        try:
                            r2_list.append(zyg_r2(u_deb,u_fin))
                        except RuntimeWarning:
                            pass
                    i_deb=i_fin+1
                    i_fin=i_deb+1
        if len(r2_list) < 2:
            # try a more exhaustive screening of SNP pairs
            r2_list=[]
            dmin=interval_list[i][0]
            dmax=interval_list[i][1]
            for chro in range(0,len(pos_list)):
                nb_snp=len(pos_list[chro])
                if nb_snp>0:
                    i_deb=0
                    i_fin=1
                    while i_fin<nb_snp:
                        while i_fin < nb_snp and (pos_list[chro][i_fin]-pos_list[chro][i_deb])<dmin:
                            i_fin+=1
                        if i_fin < nb_snp and (pos_list[chro][i_fin]-pos_list[chro][i_deb])<=dmax:
                            # compute r2
                            u_deb=geno_list[chro][:,i_deb]
                            u_fin=geno_list[chro][:,i_fin]
                            try:
                                r2_list.append(zyg_r2(u_deb,u_fin))
                            except RuntimeWarning:
                                pass
                        i_deb+=1
                        i_fin=i_deb+1
        if len(r2_list) >= 2:
            moy[i]=np.mean(np.array(r2_list,dtype='float'))
            var[i]=np.std(np.array(r2_list,dtype='float'))
    return moy,var

def distrib_ibs(pos_list,distance_list,dmax=np.inf):
    '''
    returns the probability of ibs exceeding a given distance for a list of distances.
    pos_list is a list of 1 dim arrays
    distance_list is a list of distances
    '''
    # builds a ibs length sample
    d=np.zeros(shape=1,dtype='int32')
    for chro in range(0,len(pos_list)):
        pos_temp=pos_list[chro]
        d_temp=pos_temp[1:]-pos_temp[:(len(pos_temp)-1)]
        d=np.concatenate((d,d_temp))
    if not dmax==np.inf:
        d=np.minimum(d,dmax*np.ones(shape=len(d),dtype='int32'))
    # computes the ecdf of this sample
    p=len(distance_list)
    cdf=-np.ones(shape=p,dtype='float')
    if len(d)>1:
        d=d[1:]
        for i in range(0,p):
            sel=(d>=distance_list[i])
            cdf[i]=sum(sel)
        cdf=cdf/len(d)
    return cdf



def ibs_quantiles_from_data(m,pos_list,data_type,data_list,prob_list,dmax=200000000,quantiles=False,moments=False):
    ''' Computes the quantiles of the ibs length distribution for a subset of m haplotypes or m diploid individuals

    WARNING: Set data_type to 1 if haplotypes, to 2 if genotypes
    if m==1 and data_type==2 this corresponds to ROH

    Arguments:
    m          int   nb of haplotypes or genotypes to randomly subsample
    pos_list   list(1 dim np.array)  positions of snps
    data_type  int                   if 1 data are haplotypes if 2 genotypes
    data_list  list(2-dim np.array)  haplotypes or genotypes
    prob_list  list(float)           vector of probabilities for which quantiles are computed
    dmax       int                   maximum length of ibs (eg. length of the segments), highly recommended to specify dmax
    quantiles  bool                  compute quantiles of ibs-length distribution?
    moments    bool                  compute moments of ibs-length distribution?

    Returns:
    np.concatenate((q,moms)) quantiles and moments 1 to 4th of ibs-length distrib

    Note:
    Author Simon, Flora
    '''


    q=np.array([])
    moms=np.array([])
    # builds a ibs length sample
    d=np.zeros(shape=1,dtype='int32')
    for chro in range(0,len(pos_list)):
        pos_temp=pos_list[chro]
        data_temp=data_list[chro]
        n=data_temp.shape[0]
        if m<n:
            # for each chromosome we randomly draw m haplotypes (or m genotypes)
            # so that we are not always using the same individuals for the computation
            # (it matters for real data)
            # and update count_temp and pos_temp
            subset=np.random.choice(n,size=m,replace=False)
            count_temp=np.sum(data_temp[subset,],axis=0)
            pos_temp=pos_temp[(count_temp>0)*(count_temp<(data_type*m))]
        if len(pos_temp)>1:
            d=np.concatenate((d,np.diff(pos_temp)))
        else:
            d=np.concatenate((d,dmax*np.ones(shape=1,dtype='int32')))
    d=np.minimum(d,dmax*np.ones(shape=len(d),dtype='int32'))
    # computes the quantiles and/or the moments of this sample
    if quantiles:
        q=sp.stats.mstats.mquantiles(d[1:],prob=prob_list,alphap=1,betap=1)
    if moments:
        moms=-np.ones(shape=4)
        moms[0]=np.mean(d[1:])
        moms[1]=np.std(d[1:])
        #mom3:skewness, mom4:kurtosis
        for m in range(3,5):
            moms[m-1]=np.mean(np.power((d[1:]-moms[0])/moms[1],m))
    return np.concatenate((q,moms))


def break_chr(pos_list,hap_list,dmax=2000000):
    '''
    breaks a list of long chromosomes into an equivalent list chromosomes with length lower than dmax
    to be used before ibs_quantiles in the case of real data sets with unequal chromosomes lengths
    '''
    pos_list_new=[]
    hap_list_new=[]
    for chro in range(0,len(pos_list)):
        print 'breaking chromosome '+str(chro)
        pos_temp=pos_list[chro]
        hap_temp=hap_list[chro]
        outlier_ind=(pos_temp>dmax)
        while np.sum(outlier_ind)>0:
            if np.prod(outlier_ind)==0:
                pos_list_new.append(pos_temp[np.logical_not(outlier_ind)])
                hap_list_new.append(hap_temp[:,np.logical_not(outlier_ind)])
            pos_temp=pos_temp[outlier_ind]-dmax
            hap_temp=hap_temp[:,outlier_ind]
            outlier_ind=(pos_temp>dmax)
        pos_list_new.append(pos_temp)
        hap_list_new.append(hap_temp)
    return pos_list_new,hap_list_new

def hap_to_geno(hap_list):
    '''
    transforms a list of haplotypes into a list of genotypes
    pairs of haplotypes are randomly sampled for each chromosome
    '''
    geno_list=[]
    for hap in hap_list:
        n=hap.shape[0]
        p=hap.shape[1]
        permut=np.random.permutation(n)
        geno=-np.ones(shape=(n/2,p),dtype='int32')
        for i in range(n/2):
            geno[i,:]=hap[permut[2*i],:]+hap[permut[2*i+1],:]
        geno_list.append(geno)
    return geno_list


def hap_to_geno_non_random(hap_list):
    '''transforms a list of haplotypes into a list of genotypes

    genotypes are built from 2 successive haplotypes (not randomly sampled) for each chromosome
    Note:
    Author Flora
    '''
    return [hap[0::2,:]+hap[1::2,:] for hap in hap_list]


# Implements extra statistics that are not PopSizeABC (Boitard et al 2016)
# Author: Flora Jay


def distrib_afibs(hap_list,pos_list,count_list):
    """
    Moments for length distributions of AF-IBS as defined by Theunert et al. 2012

    Arguments:
    hap_list        list(np.array[Nhap,Nsnp_seg])    haplotype data for each segment
    pos_list        list(np.array[Nsnp_seg])         positions of SNP for rach segment
    count_list      list(np.array[Nsnp_seg])         number of derived alleles at each position for each segment

    Return:
    mean_sd_afibs   np.array(mean_2,sd_2, mean_3,sd_3, ...)        (mean,sd) of afibs lengths for each category of derived alleles number 2..n
    """
    Nhap=hap_list[0].shape[0]
    afibs=[[] for der in xrange(Nhap)]
    for chro in range(0,len(pos_list)):
        afibs=afibs_fast(hap_list[chro]==1,pos_list[chro],count_list[chro],afibs)   #hap_list[chro]==1 because afibs_fast takes a boolean array as argument, not int

    mean_sd_afibs=np.zeros(shape=(len(afibs)-2)*2)
    # we don't compute afibs values for singletons (does not make sense) nor for fixed derived (because we don't simulated fixed derived, and the ones appearing because of errors added afterwards are pruned)
    i=0
    for der in xrange(2,len(afibs)):
        if len(afibs[der])>0:
            mean_sd_afibs[i]=np.mean(afibs[der])
            mean_sd_afibs[i+1]=np.std(afibs[der])
        i+=2
    return mean_sd_afibs
#np.array(list(flatten( [ (np.mean(afibs[der]),np.std(afibs[der])) if len(afibs[der])>0 else (0,0) for der in xrange(2,len(afibs)) ] )))

'''
def afibs(hap_list,pos_list,count_list):
    afibs=[[] for der in xrange(Nhap+1)]
    res=[afibs_fast(hap_list[chro]==1,pos_list[chro],count_list[chro],[[] for der in xrange(Nhap+1)])  for chro in range(0,len(pos_list))]
    return res
'''

def afibs_fast(data,posOnChrom,counts,afibs):
    """
    Compute AF-IBS as defined by Theunert et al. 2012

    Arguments:
    data          bool np.array    Rows: individuals (2 lines per individual if diploid) ; columns = snps ; True if derived
    posOnChrom    int array      positions (bp) of each polymorphism relative to its chromosome
    afibs  list(list)         list containing for each frequency f the list of lengths of afibs tracts around derived alleles at freq. f (eg afibs[2] contains list of len for sites with 2 derived alleles) for segments already analysed

    Return:
    afibs    the list is updated with the current segment's statistics

    Notes:
    An alternative version based on PBWT algorithm can be found on the gitlab, see AFIBS_BWT folder
    """
    '''
    execution time was ~ 2 s for 2Mb segment 20 haplotypes
    faster than a regular search
    it has not been tested for a larger number of individuals
    (For large number of haplotypes, the number of observed configurations might increase, and this algorithm becomes less efficient)
    An alternative version based on PBWT algorithm can be found on the gitlab, see AFIBS_BWT folder
    '''
    Nhap,Nsnp=data.shape
    ibs=np.zeros(Nsnp)
    curr_borders=dict()

    for snp in xrange(Nsnp):
        Nderived = counts[snp]
        if Nderived<=1: continue
        # Calculates a "code" corresponding to the "configuration" of each snp:
        # For each snp, the vector of 0s and 1s (ancestral and derived alleles stored for each haplo) is converted to a number in base 10
        # conf[snp] = sum_i (haplo_i * 2^i)   , (haplo_i=0 if ancestral, 1 if derived)
        # This avoid to recalculate the borders if there is a snp close by with the same configuration and for which we already calculated the borders
        conf=apply(lambda x:2**x,np.where(data[:,snp])).sum()

        if conf in curr_borders:
            # if the last calculated right border for the same conf is further right,
            # then the snp is in the same AF-IBS segment than the previous snp having this conf
            # so the borders of the segment do not change

            if curr_borders[conf][1]>snp:
                l,r = curr_borders[conf]
                if l>=0 and r!=Nsnp: afibs[Nderived].append(posOnChrom[r]-posOnChrom[l])
                continue
            else:
                # Left border cannot be before the previous left border for the same conf
                minPotentialLeft= curr_borders[conf][1]
                mask=data[:,snp]
                vec=data[mask,minPotentialLeft:snp].sum(axis=0) % Nderived
                l=np.where(vec!=0)[0][-1] + minPotentialLeft
                r,foundr=snp,False
                while not foundr and r < Nsnp-1 :
                    r+=1
                    foundr = data[mask,r].sum() % Nderived
                if not foundr:
                    # No right border in analysed segment
                    curr_borders[conf]=[l,Nsnp]
                else:
                    curr_borders[conf]=[l,r]
                    afibs[Nderived].append(posOnChrom[r]-posOnChrom[l])
                continue

        # NOT in CURR_BORDERS
        mask=data[:,snp]
        if Nderived==2:
        # It is faster when Nderived is small to compute vec
            vec=data[mask,:].sum(axis=0) % Nderived
            try:
                l=np.where(vec[:snp]!=0)[0][-1]
            except IndexError:
                l=-1
            try:
                r=np.where(vec[snp+1:]!=0)[0][0]+snp+1
            except IndexError:
                r=Nsnp
            curr_borders[conf]=[l,r]
            if l>=0 and r!=Nsnp:
                afibs[Nderived].append(posOnChrom[r]-posOnChrom[l])
            continue

        l,foundl=snp,False
        while not foundl and l > 0:
            l-=1
            foundl = data[mask,l].sum() % Nderived
        l -= int(not foundl)

        r,foundr=snp,False
        while not foundr and r < Nsnp-1 :
            r+=1
            foundr = data[mask,r].sum() % Nderived
        r+=int(not foundr)

        curr_borders[conf]=[l,r]
        if l>=0 and r!=Nsnp:
            afibs[Nderived].append(posOnChrom[r]-posOnChrom[l])

    return afibs



def classical_stats(nhaplo,count_list):
    """
    For each segment: computes heterozygosity per site, diversity (pairwise differences) per site, and Tajima's D for the segment

    Arguments:
    nhaplo          int                              total number of haplotypes
    count_list      list(np.array[Nsnp_seg])         number of derived alleles at each position for each segment

    Return:
    mean and std of these statistics
    """

    Nhap=np.float(nhaplo)
    a1 = np.sum([1.0/i for i in xrange(1,nhaplo)])
    a2 = np.sum([1.0/i**2 for i in xrange(1,nhaplo)])
    c1 = (Nhap+1)/(3*(Nhap-1)) - 1/a1   #b1-1/a1
    c2 = 2*(Nhap**2+Nhap+3)/(9*Nhap*(Nhap-1)) - (Nhap+2)/(a1*Nhap) + a2/a1**2

    all_H,all_PI,all_D = [],[],[]
    for counts in count_list:
        Nsnp=np.float(counts.shape[0])
        # Expected heterozygosity (at each site) for snp data Arlequin 8.1.1.2 p.115
        all_H.extend( 2.0/(Nhap-1) * (counts-counts**2 / Nhap)  )

        # Mean number of pariwise difference  (at each site) for snp data Arlequin 8.1.2.1 p.116
        PI =  2.0 / (Nhap * (Nhap-1)) * (counts * (Nhap-counts))
        all_PI.extend(PI)
        theta_pi = sum(PI)
        # Other estimate of theta :
        theta_s=Nsnp/a1
        #var_theta_s = (a1**2 * Nsnp + a2 * Nsnp**2) / (a1**2 * (a1**2 + a2) )
        #var_PI= (3*Nhap*(Nhap+1)*PI + 2*(Nhap**2+Nhap+3)*PI**2) / (11*(Nhap**2-7*Nhap+6))
        #var_theta_pi= (3*Nhap*(Nhap+1)*theta_pi + 2*(Nhap**2+Nhap+3)*theta_pi**2) / (11*(Nhap**2-7*Nhap+6))

        # Tajima D, formula from Tajim's paper (1989)
        all_D.append( (theta_pi - theta_s) / np.sqrt(c1/a1 *Nsnp + (c2/(a1**2+a2)) *Nsnp*(Nsnp-1)) )

    return np.array([np.mean(all_H),np.std(all_H),np.mean(all_PI),np.std(all_PI),np.mean(all_D),np.std(all_D)])



def het_one_win(dataslice,Nhap):
    '''
    Compute haplotypic heterozygosity of a given window

    Arguments:
    dataslice      np.array       subset of the data corresponding to a given window of the sequence
    Nhap           int            total number of haplotypes

    Return:
    het            float          haplotypic heterozygosity of dataslice
    '''

    haplos=[''.join([`num` for num in dataslice[i,:]]) for i in xrange(Nhap)]
    tab=Counter(haplos)
    return 1.0-sum([x**2 for x in tab.values()])/float(Nhap)**2

def haplo_win(hap_list,pos_list,win_size,L=2000000):
    '''
    Compute haplotypic heterozygosity in windows sliding aloong the genome and return mean and variance

    Arguments:
    hap_list        list(np.array[Nhap,Nsnp_seg])    haplotype data for each segment
    pos_list        list(np.array[Nsnp_seg])         positions of SNP for rach segment
    win_size        int                              lengh of the sliding windows considered as haplotypes (bp)
    L               int                              length of each simulated segment (bp)

    Return:
    mean, std       float,float                      mean and standard deviation of haplotypic heterozygosity
    '''


    Nhap=hap_list[0].shape[0]
    L=int(L)
    win_size=int(win_size)
    hetsall=[]
    for chro in xrange(len(pos_list)):
        chunks=[bisect.bisect(pos_list[chro],x) for x in range(0,L,win_size)]
        #print "chunks", chunks
        #print "pos_list chro %d "%chro, pos_list[chro]
        hets=[het_one_win(hap_list[chro][:,chunks[i]:chunks[i+1]], Nhap) for i in xrange(len(chunks)-1)]
        #print "hets chr%d: "%chro, hets
        hetsall.extend(hets)
        #tot_nb_snp+=len(pos_list[chro])
    Nhap=np.float(Nhap)
    return np.array((Nhap/(Nhap-1.0) * np.mean(hets),  Nhap/(Nhap-1.0) * np.std(hets)))



def het_one_win_durbin(x):
    '''
    Based on the algorithm described by Durbin in Efficient haplotype matching and storage using the
    Positional Burrows-Wheeler Transform (PBWT). Bioinformatics 2014
    This is giving correct results (compared to het_one_win(...) ) but does not decrease computational time
    (but likely depends on the number of haplotypes analyzed)
    For this reason I am not currently using it
    '''
    Nhap,Nsnp=x.shape
    acurr= range(Nhap)
    dcurr= [0]*Nhap
    for k in xrange(Nsnp):
        p,q=k+1,k+1
        a,b,d,e = [[] for i in xrange(4)]
        for i in xrange(Nhap):
            if dcurr[i]>p:
                p=dcurr[i]
            if dcurr[i] >q:
                q=dcurr[i]
            if x[acurr[i],k]==0:
                a.append(acurr[i])
                d.append(p)
                p=0
            else:
                b.append(acurr[i])
                e.append(q)
                q=0
        acurr=a+b
        dcurr=d+e
    haplofreq=[]
    for d in dcurr:
        if d==0:
            haplofreq[-1]+=1
        else:
            haplofreq.append(1)

    return 1.0-sum([freq**2 for freq in haplofreq])/float(Nhap)**2


def haplo_win_durbin(hap_list,pos_list,win_size,L=2000000):
    Nhap=hap_list[0].shape[0]
    L=int(L)
    win_size=int(win_size)
    hetsall=[]
    for chro in xrange(len(pos_list)):
        chunks=[bisect.bisect(pos_list[chro],x) for x in range(0,L,win_size)]
        hets=[het_one_win_durbin(hap_list[chro][:,chunks[i]:chunks[i+1]]) for i in xrange(len(chunks)-1)]
        #tot_nb_snp+=len(pos_list[chro])
        hetsall.extend(hets)
    # only last chrom return Nhap/(Nhap-1.0) * np.mean(hets),  Nhap/(Nhap-1.0) * np.std(hets)
    return Nhap/(Nhap-1.0) * np.mean(hetsall),  Nhap/(Nhap-1.0) * np.std(hetsall)




def spatial_histo_fast_unpol(polarized_list,pos_list,count_list,M,dmax=np.inf):
    '''
    Fast version of spatial_histo
    This is the correct implementation of dist
    Author Flora

    Arguments:
    polarized_list  list(np.array[Nsnp_seg])         polarized or not? for SNP for each segment
    pos_list        list(np.array[Nsnp_seg])         positions of SNP for each segment
    count_list      list(np.array[Nsnp_seg])         number of derived alleles at each position for each segment

    Returns :
    - the site frequency spectrum from 1 to M (percentages, sum to 1)
    - the variance of the distance between sites with count i, for all i. This variance needs to be multiplied by the overall proportion of SNPs.
    positions of vector pos are assumed to be SORTED within each chromosome
    '''
    histo=np.zeros(shape=M,dtype='float')
    nb_snp=0
    d = [[] for i in range(1,M+1) ]
    for chro in range(0,len(pos_list)):
        pos= [[] for i in range(1,M+1) ]
        for snp in xrange(pos_list[chro].shape[0]):
            if not polarized_list[chro][snp]:
                continue
            i=count_list[chro][snp]
            try:
                histo[i-1]+=1
                pos[i-1].append(pos_list[chro][snp])
            except IndexError:
                continue

        [d[i-1].append(x) for i in range (1,M+1) for x in np.diff(pos[i-1]) if x<=dmax]

    # for each frequency, compute the std of the distance list, after removing distance longer than dmax
    # if no distances deviation set to 0 (was -1 before)
    dist=np.asarray([np.std(d_at_freq) if len(d_at_freq)>1 else 0.0 for d_at_freq in d])

    return histo/np.sum(histo),dist




def distrib_afibs_unpol(polarized_list,hap_list,pos_list,count_list):
    """
    Moments for length distributions of AF-IBS as defined by Theunert et al. 2012

    Arguments:
    hap_list        list(np.array[Nhap,Nsnp_seg])    haplotype data for each segment
    pos_list        list(np.array[Nsnp_seg])         positions of SNP for rach segment
    count_list      list(np.array[Nsnp_seg])         number of derived alleles at each position for each segment

    Return:
    mean_sd_afibs   np.array(mean_2,sd_2, mean_3,sd_3, ...)        (mean,sd) of afibs lengths for each category of derived alleles number 2..n
    """
    Nhap=hap_list[0].shape[0]
    afibs=[[] for der in xrange(Nhap)]
    for chro in range(0,len(pos_list)):
        afibs=afibs_fast_unpol(polarized_list[chro],hap_list[chro]==1,pos_list[chro],count_list[chro],afibs)   #hap_list[chro]==1 because afibs_fast takes a boolean array as argument, not int

    mean_sd_afibs=np.zeros(shape=(len(afibs)-2)*2)
    # we don't compute afibs values for singletons (does not make sense) nor for fixed derived (because we don't simulated fixed derived, and the ones appearing because of errors added afterwards are pruned)
    i=0
    for der in xrange(2,len(afibs)):
        if len(afibs[der])>0:
            mean_sd_afibs[i]=np.mean(afibs[der])
            mean_sd_afibs[i+1]=np.std(afibs[der])
        i+=2
    return mean_sd_afibs


def afibs_fast_unpol(is_polarized,data,posOnChrom,counts,afibs):
    """
    Compute AF-IBS as defined by Theunert et al. 2012
    Some SNPs might be unpolorized (ie we do not not which allele is the derived one)

    Arguments:
    is_polarized  bool np.array(Nsnp)  'polarization status' of each SNP
    data          bool np.array    Rows: individuals (2 lines per individual if diploid) ; columns = snps ; True if derived
    posOnChrom    int array      positions (bp) of each polymorphism relative to its chromosome
    afibs  list(list)         list containing for each frequency f the list of lengths of afibs tracts around derived alleles at freq. f (eg afibs[2] contains list of len for sites with 2 derived alleles) for segments already analysed

    Return
    afibs    the list is updated with the current segment's statistics

    """
    '''
    execution time was ~ 2 s for 2Mb segment 20 haplotypes
    faster than a regular search
    it has not been tested for a larger number of individuals
    (For large number of haplotypes, the number of observed configurations might increase, and this algorithm becomes less efficient)
    '''
    Nhap,Nsnp=data.shape
    ibs=np.zeros(Nsnp)
    curr_borders=dict()

    for snp in xrange(Nsnp):
        if not is_polarized[snp]: continue
        Nderived = counts[snp]
        if Nderived<=1: continue
        # Calculates a "code" corresponding to the "configuration" of each snp:
        # For each snp, the vector of 0s and 1s (ancestral and derived alleles stored for each haplo) is converted to a number in base 10
        # conf[snp] = sum_i (haplo_i * 2^i)   , (haplo_i=0 if ancestral, 1 if derived)
        # This avoid to recalculate the borders if there is a snp close by with the same configuration and for which we already calculated the borders
        conf=apply(lambda x:2**x,np.where(data[:,snp])).sum()

        if conf in curr_borders:

            if curr_borders[conf][1]>snp:
                # if the last calculated right border for the same conf is further right,
                # then the snp is in the same AF-IBS segment than the previous snp having this conf
                # so the borders of the segment do not change
                l,r = curr_borders[conf]
                if l>=0 and r!=Nsnp: afibs[Nderived].append(posOnChrom[r]-posOnChrom[l])
                continue
            else:
                # Need to look for AF-IBS segment border
                # Left border cannot be before the previous left border for the same conf
                minPotentialLeft= curr_borders[conf][1]
                mask=data[:,snp]
                vec=data[mask,minPotentialLeft:snp].sum(axis=0) % Nderived
                l=np.where(vec!=0)[0][-1] + minPotentialLeft
                r,foundr=snp,False
                while not foundr and r < Nsnp-1 :
                    r+=1
                    foundr = data[mask,r].sum() % Nderived
                if not foundr:
                    # No right border in analysed segment
                    curr_borders[conf]=[l,Nsnp]
                else:
                    curr_borders[conf]=[l,r]
                    afibs[Nderived].append(posOnChrom[r]-posOnChrom[l])
                continue

        # NOT in CURR_BORDERS
        mask=data[:,snp] #mask=boolean ?sequence carries the derived allele
        if Nderived==2:
            # It is faster when Nderived is small to compute vec
            # vec[i]=0 if the site i is fixed (derived or ancestral) for  the subset of sequences defined by mask
            # we are looking for the next left and next right "1" in vec to define the borders around the site "snp"
            vec=data[mask,:].sum(axis=0) % Nderived
            try:
                l=np.where(vec[:snp]!=0)[0][-1]
            except IndexError:
                l=-1
            try:
                r=np.where(vec[snp+1:]!=0)[0][0]+snp+1
            except IndexError:
                r=Nsnp
            curr_borders[conf]=[l,r]
            if l>=0 and r!=Nsnp:
                afibs[Nderived].append(posOnChrom[r]-posOnChrom[l])
            continue

        # CASE wher Nderived > 2
        # we do not compute vec for all sites but look for the border by moving away from the snp 1 site by 1
        l,foundl=snp,False
        while not foundl and l > 0:
            l-=1
            foundl = data[mask,l].sum() % Nderived
        l -= int(not foundl)

        r,foundr=snp,False
        while not foundr and r < Nsnp-1 :
            r+=1
            foundr = data[mask,r].sum() % Nderived
        r+=int(not foundr)

        curr_borders[conf]=[l,r]
        if l>=0 and r!=Nsnp:
            afibs[Nderived].append(posOnChrom[r]-posOnChrom[l])

    return afibs
