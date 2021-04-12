#!/usr/bin/env python
# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# In[7]:


from copy import deepcopy
import datetime, itertools
import pandas as pd
import numpy as np

from sklearn.linear_model import LinearRegression

class Model():
    def __init__(self, factors, coefs, intercept):
        self.factors = factors
        self.model = LinearRegression()
        self.model.coef_ = np.array(coefs)
        self.model.intercept_ = intercept
    
    def predict(self, X):
        return self.model.predict(X)

mod = Model(['workplaces', 'transit', 'retail'],
            [-0.03860186,  0.16343141, 0.08733047],
            22.35632122147008
           )


pd.set_option('display.max_columns', 500)


# In[19]:


mod.model.coef_[mod.factors.index('transit')]


# In[3]:


# df = pd.read_csv('Global_Mobility_Report.csv') # old file just for testing
# df


# In[28]:


# LONG

print("Downloading data from Google..."
import gzip
import urllib.request

request = urllib.request.Request(
    "https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv",
    headers={
        "Accept":"text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
        "Accept-Encoding":"gzip, deflate, br",
        "Accept-Language":"en-US,en;q=0.5",
        "Cache-Control":"no-cache",
        "Connection":"keep-alive",
        "DNT":"1",
        "Host":"www.gstatic.com",
        "Pragma":"no-cache",
        "TE":"Trailers",
        "Upgrade-Insecure-Requests":"1",
        "User-Agent":"Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:82.0) Gecko/20100101 Firefox/82.0",
    }
)
response = urllib.request.urlopen(request)
gzipFile = gzip.GzipFile(fileobj=response)


# In[29]:


europe_country_codes = ['AL / Albania',
                        'AD / Andorra',
                        'AM / Armenia',
                        'AT / Austria',
                        'BY / Belarus',
                        'BE / Belgium',
                        'BA / Bosnia and Herzegovina',
                        'BG / Bulgaria',
                        'CH / Switzerland',
                        'CY / Cyprus',
                        'CZ / Czech Republic',
                        'DE / Germany',
                        'DK / Denmark',
                        'EE / Estonia',
                        'ES / Spain',
                        'FO / Faeroe Islands',
                        'FI / Finland',
                        'FR / France',
                        'GB / United Kingdom',
                        'GE / Georgia',
                        'GI / Gibraltar',
                        'GR / Greece',
                        'HU / Hungary',
                        'HR / Croatia',
                        'IE / Ireland',
                        'IS / Iceland',
                        'IT / Italy',
                        'LT / Lithuania',
                        'LI / Liechtenstein',
                        'LU / Luxembourg',
                        'LV / Latvia',
                        'MC / Monaco',
                        'MK / Macedonia',
                        'MT / Malta',
                        'NO / Norway',
                        'NL / Netherlands',
                        'PO / Poland',
                        'PL / Poland',
                        'PT / Portugal',
                        'RO / Romania',
                        'RU / Russian Federation',
                        'SE / Sweden',
                        'SI / Slovenia',
                        'SK / Slovakia',
                        'SM / San Marino',
                        'TR / Turkey',
                        'UA / Ukraine',
                        'VA / Vatican City State'
                       ]
europe_country_codes = [cc.split(' / ')[0] for cc in europe_country_codes]

def filter_europe(df):
    df.date = pd.to_datetime(df.date)
    res = df[(df.country_region_code.isin(europe_country_codes)) & (df.sub_region_1.isnull())]
    res = res.drop(columns=['country_region',
                            'sub_region_1',
                            'sub_region_2',
                            'metro_area',
                            'iso_3166_2_code',
                            'census_fips_code',
                            'place_id',
                           ]).set_index('date').sort_index()
    return res.rename(columns=lambda c:c.split("_")[0])


# In[8]:


eu = filter_europe(pd.read_csv(gzipFile))
print("File from Google downloaded")
eu


# In[9]:


def get_country_means(eu, country_code):
    c = eu[eu.country == country_code]
    if len(c) == 0:
        return None
    day = c.index[-1].dayofweek # dayofweek je číslo dne, pondělí = 0, neděle = 6
    if day != 6:
        c = c.iloc[:-day-1]
    return c.resample('W').mean()


# In[10]:


cz = get_country_means(eu, 'CZ')
cz


# In[30]:


def predict_country(c, mod):
    tmp = c[mod.factors].dropna()
    index = tmp.index
    X = tmp[mod.factors].values
    y = mod.model.predict(X)
    c['contacts'] = pd.Series(y, index)
    return c

def predict_li(c, mod):
    tmp = c['transit'].dropna()
    index = tmp.index
    X = tmp.values
    y = X * mod.model.coef_[mod.factors.index('transit')] + mod.model.intercept_
    c['contacts'] = pd.Series(y, index)
    return c


# In[31]:


predict_country(cz, mod)


# In[34]:


for cc in europe_country_codes:
    c = get_country_means(eu, cc)
    if c is not None:
        print("Predicting contacts for %s" % cc)
        if cc == "LI":
            predict_li(c, mod).to_csv('country_data/%s.csv' % cc)
        else:
            predict_country(c, mod).to_csv('country_data/%s.csv' % cc)


# In[ ]:




