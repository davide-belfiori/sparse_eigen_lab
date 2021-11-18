# Importazione dei file necessari
load thermomech_TC.mat;
arnoldi
# Generazione del vettore b di partenza
b = randn(size(Problem.A)(2), 1);
# Esecuzione del metodo di Arnoldi
n = 10
[Q, H] = arnoldi_iter(Problem.A, b, n);
# controllo ortogonalità
check_orth = norm(eye(n+1) - Q'*Q)
# controllo errore residuo
residual = norm(Problem.A*Q(:,1:n)-Q*H)

H_n = H(1:n, :)

# calcolo dei valori di Ritz
n_val = 10
ritz = sort(eig(H_n), "descend")(1:n_val)

# calcolo digli autovalori 'reali'
eg_max = sort(eigs(Problem.A, n_val, "lm"), "descend")

# errore di approssimazione
approx_err = abs(eg_max - ritz)

# Variazione di n

n = 20
[Q, H] = arnoldi_iter(Problem.A, b, n);
check_orth = norm(eye(n+1) - Q'*Q)
residual = norm(Problem.A*Q(:,1:n)-Q*H)
H_n = H(1:n, :);
ritz_1 = sort(eig(H_n), "descend")(1:n_val)
approx_err_1 = abs(eg_max - ritz_1)

n = 50
[Q, H] = arnoldi_iter(Problem.A, b, n);
check_orth = norm(eye(n+1) - Q'*Q)
residual = norm(Problem.A*Q(:,1:n)-Q*H)
H_n = H(1:n, :);
ritz_2 = sort(eig(H_n), "descend")(1:n_val)
approx_err_2 = abs(eg_max - ritz_2)

# plot dei risultati
#plot(eg_max,'rx',ritz,'bo',ritz_1,'ro',ritz_2,'ko')
#t = title("Distribuzione");
#h = legend("autoval.", "n = 5","n = 10", "n = 50");
#legend (h, "location", "northeastoutside");
#set (gca, "fontsize", 20);


plot(approx_err,'b',approx_err_1,'r',approx_err_2,'k')
t = title("Errore di approssimazione (pos. def.)");
h = legend("n = 10","n = 20", "n = 50");
legend (h, "location", "northeastoutside");
set (gca, "fontsize", 20);