import numpy as np
import pandas as pd
import time as t
import random

#importing homemade modules
import predictor as p
import liveOperator as lo
import orderer as o

#USD value of funds actively managed by the bot
principalValue = 10000

#We will wait a certain number of milliseconds after booting up the bot
delayStart = int((10/60)*3600000) #10 minute delay before starting up again

#Below are 5 cryptos with high trading volume.  These may or may not be cryptos that I am actually trading!
tradingInstruments = ["BTC","ETH","USDT","XRP","BNB"]

startTime = int(t.time()*1000) + delayStart
prevDH=None
lastHour = None

#The bot reassesses its holdings every hour, forever
while True:
    rightNow = int(t.time()*1000) #ms UNIX time
    hour = int(rightNow/3600000)
    delayOver = (rightNow >= startTime)
    notStartedYet = (lastHour is None)
    if (delayOver and (notStartedYet or ((rightNow - lastHour*3600000) >= (3600000 + r)))): #if it's been more than an hour since last time
        try:
            lastHour = hour
            """
            r is chosen as a random number here
            """

            #Predictor object streams live market data and predicts future price changes
            pea = p.Predictor()
            pea.appendPrediction()
        
            #Orderer object handles all communications with the exchange (such as buy/sell orders)
            orwell = o.Orderer()

            if prevDH is None:
                prevDH = pd.Series(orwell.getBalance(tradingInstruments))
        
            #"live operator" object decides actions to take given predictions
            leo = lo.LiveOperator(pea, orwell, tradingInstruments, principalValue)
            if prevDH is not None:
                leo.loadPreviousDesiredHoldings(prevDH)
            leo.submitOrders()
            leo.appendDesiredHoldings()
            prevDH = leo.desiredHoldings
        except Exception as e: #if there is any error, we don't want it to stop the whole program
            print("Uncaught error! "+str(rightNow))
            print(e)
            #For what it's worth, this exception has never been triggered before!
    t.sleep(15) #sleep for 15 seconds

