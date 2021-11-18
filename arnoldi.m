1;
function [Q, H] = arnoldi_iter(A, b, n, eps = 10e-16, reort = 1)
  # Calcola una base ortonormale dell'(n+1)-esimo sottospazio di Krylov 
  # formanto dalla matrice A e dal vettore b.
  
  % Inizializzazione delle matrici H e Q
  H = zeros(n + 1, n);
  Q = zeros(size(A)(1), n+1);
  
  % Inizializzazione della prima colonna di Q
  Q(:,1) = b / norm(b); # primo vettore dello spazio di Krylov normalizzato
  
  for j = 1 : n
    % Generazione di un nuovo vettore della base
    v = A * Q(:,j);
    % Ortogonalizzazione
    for k = 1 : j
      % Sottrazione delle proiezioni sui vettori precedenti
      H(k, j) = v' * Q(:,k);
      v = v - H(k,j) * Q(:,k);
    endfor
    if reort
       # Ri-ortogonalizzazione
       for k = 1:j          
         v = v - (v' * Q(:,k)) * Q(:,k);
       endfor
     endif
    % Normalizzazione
    H(j+1, j) = norm(v);
    if H(j+1,j) > eps
      Q(:,j+1) = v / H(j+1,j);
    else
      'Breakdown'
      return
    endif  
  endfor
endfunction


function [Q, H] = inverse_arnoldi_iter(A, b, n, eps = 10e-16,gmres_iter = 50, reort = 1)
  # Calcola una base ortonormale dell'(n+1)-esimo sottospazio di Krylov 
  # formanto dalla matrice A^(-1) e dal vettore b.
  
  % Inizializzazione delle matrici H e Q
  H = zeros(n + 1, n);
  Q = zeros(size(A)(1), n+1);
  
  % Inizializzazione della prima colonna di Q
  Q(:,1) = b / norm(b); # primo vettore dello spazio di Krylov normalizzato
  
  for j = 1 : n
    % Generazione di un nuovo vettore della base
    v = A \ Q(:,j);
    % Ortogonalizzazione
    for k = 1 : j
      % Sottrazione delle proiezioni sui vettori precedenti
      H(k, j) = v' * Q(:,k);
      v = v - H(k,j) * Q(:,k);
    endfor
    # Ri-ortogonalizzazione
    if reort
       for k = 1:j          
         v = v - (v' * Q(:,k)) * Q(:,k);
       endfor
    endif
    % Normalizzazione
    H(j+1, j) = norm(v);
    if H(j+1,j) > eps
      Q(:,j+1) = v / H(j+1,j);
    else
      'Breakdown'
      return
    endif  
  endfor
endfunction
