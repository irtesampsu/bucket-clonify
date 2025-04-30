from sklearn.metrics.cluster import adjusted_rand_score, adjusted_mutual_info_score

def calculate_adjusted_rand_score(labels_true: list, labels_pred: list) -> float:
    assert len(labels_true) == len(labels_pred), "[assert] labels_true and labels_pred must have same number of elements!"

    return adjusted_rand_score(labels_true=labels_true, labels_pred=labels_pred)

def calculate_adjusted_mutual_info_score(labels_true: list, labels_pred: list) -> float:
    assert len(labels_true) == len(labels_pred), "[assert] labels_true and labels_pred must have same number of elements!"

    return adjusted_mutual_info_score(labels_true=labels_true, labels_pred=labels_pred)

if __name__ == "__main__":
    labels_true, labels_pred = [1, 1, 0, 1, 0], [0, 0, 1, 0, 1]

    print(f"Calculated Adjusted Rand Index (ARI) = {calculate_adjusted_rand_score(labels_true=labels_true, labels_pred=labels_pred)}")
    print(f"Calculated Adjusted Mutual Information (AMI) = {calculate_adjusted_mutual_info_score(labels_true=labels_true, labels_pred=labels_pred)}")