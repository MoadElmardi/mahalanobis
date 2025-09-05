use std::vec;
use tfhe::{prelude::*, ClientKey, FheInt32Id, PublicKey, Tag};
use tfhe::{ConfigBuilder, generate_keys, set_server_key, FheInt32, CpuFheInt32Array};
use tfhe::integer::ciphertext::BaseSignedRadixCiphertext;
use tfhe::shortint::Ciphertext;

fn mean_vector(d: CpuFheInt32Array, n_cols: usize, public_key: PublicKey) -> CpuFheInt32Array {
    // Nombre de lignes (d.shape()[0] renvoit le nombre d'éléments de d)
    let n_rows = d.shape()[0] / n_cols;
    
    // Calcul des sommes des éléments des colonnes
    let mut sum_array: CpuFheInt32Array = d.slice(&[0..n_cols]).into_owned();

    for i in 1..n_rows {
        let slice = d.slice(&[i*n_cols..(i+1)*n_cols]);
        sum_array = &sum_array + &slice;
    }

    // Création d'un tableau contenant des éléments égaux au nombre de lignes et de taille égale au nombre de colonnes
    let n_array = vec![n_rows as i32; n_cols];
    let n_encrypted_vec: Vec<FheInt32> = n_array.iter().map(|&x| FheInt32::encrypt(x, &public_key)).collect();
    let shape = vec![n_cols];
    let ct_vec: Vec<BaseSignedRadixCiphertext<Ciphertext>> = n_encrypted_vec
    .into_iter()
    .map(|x| x.into_raw_parts().0)
    .collect();
    let n_array_encrypted = CpuFheInt32Array::new(ct_vec, shape);
    
    // Calcul du vecteur moyen
    let mean_array = &sum_array / &n_array_encrypted;
    return mean_array;
}

fn covariance_matrix(d: CpuFheInt32Array, n_cols: usize, public_key: PublicKey) -> CpuFheInt32Array {
    // Nombre de lignes
    let n_rows = d.shape()[0] / n_cols;

    // Vecteur moyen
    let mean = mean_vector(d.clone(), n_cols, public_key.clone());

    // Création d'une matrice chiffrée remplie par des 0
    let covariance_clair = vec![0i32; n_cols * n_cols];
    let covariance_chiffree: Vec<FheInt32> = covariance_clair.iter().map(|&x| FheInt32::encrypt(x, &public_key)).collect();
    let shape = vec![n_cols * n_cols];
    let ct_vec: Vec<BaseSignedRadixCiphertext<Ciphertext>> = covariance_chiffree
        .into_iter()
        .map(|x| x.into_raw_parts().0)
        .collect();
    let mut covariance = CpuFheInt32Array::new(ct_vec, shape);

    for i in 0..n_rows {
        // Trancher la ligne i
        let row_slice = d.slice(&[i*n_cols..(i+1)*n_cols]);
        // Calculer la différence entre la ligne i et le vecteur moyen
        let diff = &row_slice - &mean;

        // Calculer le produit extérieur entre diff et la transposée de diff
        for j in 0..n_cols {
            for k in 0..n_cols {
                let outer_product = diff.slice(&[j..j+1]) * diff.slice(&[k..k+1]);
                
                // Créer une matrice de taille[n_cols, n_cols] remplie de zéros
                let mut update_vec = vec![FheInt32::encrypt(0i32, &public_key); n_cols * n_cols];

                // Placer le outer_product à la position (j, k)
                let idx = j * n_cols + k;
                let raw_cipher = outer_product.container()[0].clone();
                let scalar = FheInt32::from_raw_parts(raw_cipher, FheInt32Id, Tag::default());
                update_vec[idx] = scalar;

                let ct_vec: Vec<BaseSignedRadixCiphertext<Ciphertext>> = update_vec
                    .into_iter()
                    .map(|x| x.into_raw_parts().0)
                    .collect();

                let update_array = CpuFheInt32Array::new(ct_vec, vec![n_cols * n_cols]);

                // Ajouter la nouvelle matrice à la matrice de covariance
                covariance = covariance + update_array;
            }
        } 
    }
    
    // Créer une matrice remplie de la valeur (n-1)
    let n_1 = n_rows - 1;
    let n_1_array = vec![n_1 as i32; n_cols * n_cols];
    let n_1_encrypted_vec: Vec<FheInt32> = n_1_array.iter().map(|&x| FheInt32::encrypt(x, &public_key)).collect();
    let n_1_shape = vec![n_cols * n_cols];
    let n_1_ct_vec: Vec<BaseSignedRadixCiphertext<Ciphertext>> = n_1_encrypted_vec
        .into_iter()
        .map(|x| x.into_raw_parts().0)
        .collect();
    let n_1_array_encrpyted = CpuFheInt32Array::new(n_1_ct_vec, n_1_shape);
    
    // Diviser la matrice de covariance par la matrice des (n-1) et obtenir le résultat correct
    let covariance = &covariance / &n_1_array_encrpyted; 
    return covariance;
}

fn gauss_jordan_inverse(mut matrix: CpuFheInt32Array, public_key: PublicKey) -> CpuFheInt32Array {
    // This function implements the Gauss-Jordan elimination method to find the inverse of the matrix

    let n = matrix.shape()[0].isqrt();
    let mut identity = vec![0i32; n * n];
    for i in 0..n {
        identity[i * n + i] = 1;   
    }

    let identity_encrypted: Vec<FheInt32> = identity.iter().map(|&x| FheInt32::encrypt(x, &public_key)).collect();
    let identity_shape = vec![n * n];
    let identity_ct_vec: Vec<BaseSignedRadixCiphertext<Ciphertext>> = identity_encrypted
        .into_iter()
        .map(|x| x.into_raw_parts().0)
        .collect();
    let mut identity_matrix = CpuFheInt32Array::new(identity_ct_vec, identity_shape);

    let zero = FheInt32::encrypt(0i32, &public_key);
    let un = FheInt32::encrypt(1i32, &public_key);

    for i in 0..n {
        // Find the pivot
        let pivot = matrix.slice(&[i*n+i..i*n+i+1]);

        if pivot.container()[0] == zero.clone().into_raw_parts().0 {
            panic!("Matrix is singular and cannot be inverted");
        }

        // Create the pivot matrix
        let pivot_value = pivot.container()[0].clone();
        let mut pivot_matrix = vec![un.clone().into_raw_parts().0; n * n];  // n = nb de colonnes
        for j in 0..n {
            pivot_matrix[i * n + j] = pivot_value.clone();
        }
        let pivot_matrix_chiffree = CpuFheInt32Array::new(pivot_matrix, vec![n * n]);

        // Divide by the pivot
        matrix = matrix / &pivot_matrix_chiffree;
        identity_matrix = identity_matrix / &pivot_matrix_chiffree;

        // Cancel out the coloumns
        let mut cancel_matrix = vec![zero.clone().into_raw_parts().0; n * n];
        let mut cancel_idendidty = vec![zero.clone().into_raw_parts().0; n * n];
        for j in 0..n {
            if j != i {
                let case = matrix.slice(&[j*n+i..j*n+i+1]);
                for k in 0..n {
                    let factor = matrix.slice(&[i*n+k..i*n+k+1]) * case.clone();
                    let factor_identity = identity_matrix.slice(&[i*n+k..i*n+k+1]) * case.clone();
                    cancel_matrix[j * n + k] = factor.container()[0].clone();
                    cancel_idendidty[j * n + k] = factor_identity.container()[0].clone();
                }
            }
        }
        let cancel_matrix_chiffree = CpuFheInt32Array::new(cancel_matrix, vec![n * n]);
        let cancel_identity_chiffree = CpuFheInt32Array::new(cancel_idendidty, vec![n * n]);

        matrix = matrix - &cancel_matrix_chiffree;
        identity_matrix = identity_matrix - &cancel_identity_chiffree;

    }

    return identity_matrix
}

fn mahalanobis_distance(d: CpuFheInt32Array, x: CpuFheInt32Array, public_key: PublicKey, client_key: ClientKey) -> CpuFheInt32Array {
    let n_cols = x.shape()[0];
    println!("n_cols: {}\n", n_cols);

    // Calculer la moyenne mu
    let mu = mean_vector(d.clone(), n_cols, public_key.clone());
    let mu_decrypted: Vec<i32> = mu.decrypt(&client_key);
    println!("Moyenne mu: {:?}\n", mu_decrypted);

    // Calculer la différence entre le point x et la moyenne mu
    let diff = &x - &mu;
    let diff_decrypted: Vec<i32> = diff.decrypt(&client_key);
    println!("Différence (x - mu): {:?}\n", diff_decrypted);

    let covariance = covariance_matrix(d, n_cols, public_key.clone());
    let covariance_decrypted: Vec<i32> = covariance.decrypt(&client_key);
    println!("Matrice de covariance: {:?}\n", covariance_decrypted);

    let covariance_inv = gauss_jordan_inverse(covariance, public_key.clone());
    let covariance_inv_decrypted: Vec<i32> = covariance_inv.decrypt(&client_key);
    println!("Matrice de covariance inverse: {:?}\n", covariance_inv_decrypted);

    // Calculer le produit (diff ^ T) * covariance_inv
    let zero = FheInt32::encrypt(0i32, &public_key);
    let mut c = CpuFheInt32Array::new(vec![zero.into_raw_parts().0; n_cols], vec![1]);
    let mut update_vec = vec![FheInt32::encrypt(0i32, &public_key); n_cols];
    for i in 0..n_cols{
        for j in 0..n_cols {
            c = c + diff.slice(&[j..j+1]) * covariance_inv.slice(&[i+n_cols*j..i+n_cols*j+1]);
        }
        let raw_cipher = c.container()[0].clone();
        let scalar = FheInt32::from_raw_parts(raw_cipher, FheInt32Id, Tag::default());
        update_vec[i] = scalar;
        c = c.clone() - c;
    }
    let ct_vec: Vec<BaseSignedRadixCiphertext<Ciphertext>> = update_vec
        .into_iter()
        .map(|x| x.into_raw_parts().0)
        .collect();
    let product = CpuFheInt32Array::new(ct_vec, vec![n_cols]);
    let product_decrypted: Vec<i32> = product.decrypt(&client_key);
    println!("Produit (diff^T * covariance_inv): {:?}\n", product_decrypted);

    // Calculer le produit scalaire (product) * diff
    for i in 0..n_cols {
        c = c + product.slice(&[i..i+1]) * diff.slice(&[i..i+1]);
    }    
    let c_decrypted: Vec<i32> = c.decrypt(&client_key);
    println!("Produit scalaire: {:?}\n", c_decrypted);

    // Calculer la distance de Mahalanobis : sqrt(diff^T * covariance_inv * diff)
    let dist = c;
    dist
}

pub fn main() {
    let config = ConfigBuilder::default().build();
    let (client_key, server_key) = generate_keys(config);
    set_server_key(server_key);
    let public_key = PublicKey::new(&client_key);

    // let x: [i32; 3] = [66, 640, 44];
    // let x: [i32; 3] = [4, 500, 40];
    // let d:Vec<i32> = vec![
    //     64, 580, 29,
    //     66, 570, 33,
    //     68, 590, 37,
    //     69, 660, 46,
    //     73, 600, 55,
    // ];
    // let d:Vec<i32> = vec![
    //     1, 100, 10,
    //     2, 300, 15,
    //     4, 200, 20,
    //     2, 600, 10,
    //     5, 100, 30,
    // ];
    // let d:Vec<i32> = vec![
    //     1, 2, 0,
    //     0, 1, 0,
    //     0, 0, 1,
    // ];
    // let d:Vec<i32> = vec![
    //     11, 50, 34,
    //     50, 1250, 205,
    //     34, 205, 110,
    // ];
    // let d:Vec<i32> = vec![
    //     1, 5, 3,
    //     5, 125, 20,
    //     3, 20, 11,
    // ];
    let d:Vec<i32> = vec![
        4, 3, 8,
        6, 2, 5,
        1, 5, 9,
    ];

    let d_encrypted = CpuFheInt32Array::try_encrypt(d.as_slice(), &client_key).unwrap();

    let covariance_inv = gauss_jordan_inverse(d_encrypted.clone(), public_key.clone());
    let covariance_inv_decrypted: Vec<i32> = covariance_inv.decrypt(&client_key);
    println!("Matrice de covariance inverse: {:?}", covariance_inv_decrypted);

    let x_encrypted = CpuFheInt32Array::try_encrypt(&x, &client_key).unwrap();

    let dist = mahalanobis_distance(d_encrypted, x_encrypted, public_key.clone(), client_key.clone());
    let dist_decrypted: Vec<i32> = dist.decrypt(&client_key);
    println!("Mahalanobis distance: {:?}", dist_decrypted);


    // let inverse = gauss_jordan_inverse(d_encrypted.clone(), public_key.clone());
    // let inverse_decrypted: Vec<i32> = inverse.decrypt(&client_key);
    // println!("Inverse matrix: {:?}", inverse_decrypted);

    // let covariance = covariance_matrix(d_encrypted, 3, public_key.clone());
    // let covariance_decrypted: Vec<i32> = covariance.decrypt(&client_key);
    // println!("Covariance matrix: {:?}", covariance_decrypted);
}