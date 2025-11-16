function SMA = compute_sma(BA, win_sma)

SMA = movmean( abs(BA(:,1))+ abs(BA(:,2)) + abs(BA(:,3)), win_sma);

end