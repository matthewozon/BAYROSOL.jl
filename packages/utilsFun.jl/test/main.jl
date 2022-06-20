    using PyPlot
using utilsFun


x = collect(-10.0:0.01:10.0);

figure(); plot(x,softMax(x)); plot(x,softMaxA(x,2.0)); plot(x,softMaxA(x,0.5)); legend(["softMax","softMax \$\\alpha=2\$","softMax \$\\alpha=0.5\$"]); title("Test softmax function")
figure(); plot(x,softMaxDeriv(x)); plot(x,softMaxDerivA(x,2.0)); plot(x,softMaxDerivA(x,0.5)); legend(["softMax","softMax \$\\alpha=2\$","softMax \$\\alpha=0.5\$"]); title("Test softmax derivative function")

y = collect(0.01:0.01:10.0)
figure(); plot(y,softMaxInv(y)); plot(y,softMaxInvA(y,2.0)); plot(y,softMaxInvA(y,0.5)); legend(["softMax\$^{-1}\$","softMax\$^{-1}\$ \$\\alpha=2\$","softMax\$^{-1}\$ \$\\alpha=0.5\$"]); title("Test softmax\$^{-1}\$ function")
figure(); semilogy(y,softMaxInvDeriv(y)); semilogy(y,softMaxInvDerivA(y,2.0)); semilogy(y,softMaxInvDerivA(y,0.5)); legend(["softMax\$^{-1}\$","softMax\$^{-1}\$ \$\\alpha=2\$","softMax\$^{-1}\$ \$\\alpha=0.5\$"]); title("Test softmax\$^{-1}\$ function")


figure(); plot(x,logistic(x)); plot(x,logistic(x,-1.0,0.0,1.0)); plot(x,logistic(x,1.0,2.0,1.0)); title("Test logistic function")
figure(); plot(x,logisticDeriv(x)); plot(x,logisticDeriv(x,-1.0,0.0,1.0)); plot(x,logisticDeriv(x,1.0,2.0,1.0)); title("Test logistic derivative function")

y1 = collect(0.01:0.01:0.99);
y2 = -1.0.+y1;
y3 = 1.0.+y1;

figure(); plot(y1,logisticInv(y1,0.0,1.0)); plot(y2,logisticInv(y2,-1.0,0.0)); plot(y3,logisticInv(y3,1.0,2.0));  title("Test logistic\$^{-1}\$ function")
