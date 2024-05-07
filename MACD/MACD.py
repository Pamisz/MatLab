import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def simulation(file, dates, prices, buy_signals, sell_signals):
    money = 0 
    product = 1000
    first_price =  product*prices[0]    
    print("\nPoczątkowy kapitał: 1000 jednostek produktu\nWartość początkowa: ", first_price, "\n")
    for i in file.index:
        if buy_signals[i] and product!= 1000:
            product = money/prices[i]
            print("Nabyto",product, " jednostek za ", money," w dniu ", dates[i])
            money = 0


        if sell_signals[i]:
            money = product*prices[i]
            print("Sprzedano",product, " jednostek za ", money, " w dniu ",dates[i])
            product = 0

    if product != 0:
        money = product*prices[i]
    return money, first_price

def calculate_and_plot_MACD(prices, file):
    short_ema = prices.ewm(span=12).mean()
    long_ema = prices.ewm(span=26).mean()
    macd = short_ema- long_ema
    signal = macd.ewm(span=9).mean()

    sell_signals = (macd.shift(1) > signal.shift(1)) & (macd < signal)
    buy_signals = (macd.shift(1) < signal.shift(1)) & (macd > signal)

    plt.plot(file.index, macd, label='MACD', color='blue')
    plt.plot(file.index, signal, label='SIGNAL', color='black')

    plt.scatter(file.index[buy_signals], macd[buy_signals], marker='.', color='green', label='Sygnał kupna')
    plt.scatter(file.index[sell_signals], macd[sell_signals], marker='.', color='red', label='Sygnał sprzedaży')

    return buy_signals, sell_signals

#Odczytanie wartości dla złota
file = pd.read_csv("gold.csv")
dates = file["Data"]
prices = file["Zamkniecie"]

#Wykres dla złota
plt.ylabel("Wartości [$/ozt]")
plt.xlabel("Indeksy akcji")
plt.title("Ceny akcji złota")
plt.plot(file.index, prices, label='Dane historyczne', color="blue")
plt.legend()
plt.savefig('akcje_gold.png', format='png')
plt.show()

#MACD
buy_signals, sell_signals = calculate_and_plot_MACD(prices,file)
plt.ylabel("Wartości [$/ozt]")
plt.xlabel("Indeksy akcji")
plt.title("Analiza MACD dla cen akcji złota")
plt.legend()
money, first_price = simulation(file, dates, prices, buy_signals, sell_signals)
print("\nFinanse po inwestycji:", money, "$\nZysk w ciągu 4 lat:", money - first_price, "$")
plt.savefig('macd_gold.png', format='png')
plt.show()
print("\n")


#Odczytanie wartości dla wig20
file = pd.read_csv("wig20.csv")
dates = file["Data"]
prices = file["Zamkniecie"]

#Wykres dla wig20
plt.ylabel("Wartości [PLN]")
plt.xlabel("Indeksy akcji")
plt.title("Ceny akcji wig20")
plt.plot(file.index, prices, label='Dane historyczne', color="blue")
plt.legend()
plt.savefig('akcje_wig20.png', format='png')
plt.show()

#MACD
buy_signals, sell_signals = calculate_and_plot_MACD(prices,file)
plt.ylabel("Wartości [PLN]]")
plt.xlabel("Indeksy akcji")
plt.title("Analiza MACD dla cen akcji wig20")
plt.legend()
money, first_price = simulation(file, dates, prices, buy_signals, sell_signals)
print("\nFinanse po inwestycji:", money, "PLN\nZysk w ciągu 4 lat:", money - first_price, "PLN")
plt.savefig('macd_wig20.png', format='png')
plt.show()