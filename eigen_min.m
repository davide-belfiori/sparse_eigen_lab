# Importazione dei file necessari
load thermal2.mat;
arnoldi
# Generazione del vettore b di partenza
b = randn(size(Problem.A)(2), 1);
# Esecuzione del metodo di Arnoldi
n = 15
[Q, H] = inverse_arnoldi_iter(Problem.A, b, n);
# controllo ortogonalità
check_orth = norm(eye(n+1) - Q'*Q)

H_n = H(1:n, :)

# calcolo dei valori di Ritz
n_val = 3
ritz = sort(eig(H_n), "ascend")(1:n_val);
ritz = flip(ritz)

# calcolo digli autovalori della matrice originale
eg_min = eigs(Problem.A, n_val, "sm")

# errore di approssimazione
approx_err = abs(eg_min - ritz)

# Variazione di n

n = 18
[Q, H] = inverse_arnoldi_iter(Problem.A, b, n);
check_orth = norm(eye(n+1) - Q'*Q)
H_n = H(1:n, :);
ritz_1 = sort(eig(H_n), "ascend")(1:n_val);
ritz_1 = flip(ritz_1)
approx_err_1 = abs(eg_min - ritz_1)

n = 20
[Q, H] = inverse_arnoldi_iter(Problem.A, b, n);
check_orth = norm(eye(n+1) - Q'*Q)
H_n = H(1:n, :);
ritz_2 = sort(eig(H_n), "ascend")(1:n_val);
ritz_2 = flip(ritz_2)
approx_err_2 = abs(eg_min - ritz_2)

# plot dei risultati
#plot(eg_min,'rx',ritz,'bo',ritz_1,'ro',ritz_2,'ko')
#t = title("Distribuzione");
#h = legend("autoval.", "n = 15","n = 18", "n = 20");
#legend (h, "location", "northeastoutside");
#set (gca, "fontsize", 20);


plot(approx_err,'b',approx_err_1,'r',approx_err_2,'k')
t = title("Errore di approssimazione");
h = legend("n = 15","n = 18", "n = 20");
legend (h, "location", "northeastoutside");
set (gca, "fontsize", 20);