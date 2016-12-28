
clear all;

WallMobGlobals;

% UF
muf_f();
muf_g();

figure;
x=MUF_f(:,1);
semilogy(x,MUF_f(:,2), x,MUF_f(:,3),x,MUF_f(:,4));
axis([0 1 10^0 10^4]);
legend('f_1^{UF}','f_3^{UF}','f_5^{UF}');

figure;
x=MUF_g(:,1);
semilogy(x,MUF_g(:,2), x,MUF_g(:,3),x,MUF_g(:,4));
axis([0 1 10^0 10^4]);
legend('g_1^{UF}','g_3^{UF}','g_5^{UF}');

% UL
mul_f();

figure;
x = MUL_f(:,1);
plot(x, MUL_f(:,2));
axis([0 1 -0.5 0.5]);
legend('f_2^{UL}');

figure;
x = MUL_f(:,1);
semilogy(x, MUL_f(:,3));
axis([0 0.5 10^-2 10^4]);
legend('f_4^{UL}');

% US
mus_f();
mus_g();

figure;
x = MUS_f(:,1);
semilogy(x, MUS_f(:,2), x, MUS_f(:,3), x, MUS_f(:,4));
axis([0 0.5 10^-1 10^4]);
legend('f_2^{US}', 'f_4^{US}', 'f_6^{US}');

figure;
x = MUS_g(:,1);
semilogy(x, MUS_g(:,2), x, MUS_g(:,3), x, MUS_g(:,4));
axis([0 0.5 10^-1 10^4]);
legend('g_2^{US}', 'g_4^{US}', 'g_6^{US}');


% OL
mol_fg();

figure;
x = MOL_fg(:,1);
semilogy(x, MOL_fg(:,2), x, MOL_fg(:,3));
axis([0 1 10^0 10^4]);
legend('f_3^{OL}', 'g_3^{OL}');


% OS
mos_f();

figure;
x = MOS_f(:,1);
semilogy(x, MOS_f(:,2), x, MOS_f(:,3));
axis([0 1 10^0 10^5]);
legend('f_3^{OS}', 'f_5^{OS}');

% ES
mes_f();
mes_g();
mes_h();

figure;
x = MES_f(:,1);
semilogy(x,MES_f(:,2), x,MES_f(:,3), x,MES_f(:,4));
axis([0 1 10^-1 10^5]);
legend('f_3^{ES}', 'f_5^{ES}', 'f_7^{ES}');

figure;
x = MES_g(:,1);
semilogy(x,MES_g(:,2), x,MES_g(:,3), x,MES_g(:,4));
axis([0 1 10^0 10^5]);
legend('g_3^{ES}', 'g_5^{ES}', 'g_7^{ES}');

figure;
x = MES_h(:,1);
semilogy(x,MES_h(:,2), x,MES_h(:,3), x,MES_h(:,4));
axis([0 1 10^0 10^5]);
legend('h_3^{ES}', 'h_5^{ES}', 'h_7^{ES}');





