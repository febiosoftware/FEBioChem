resource "github_actions_secret" "aws-access-key-id" {
  repository      = var.repo_name
  secret_name     = "AWS_ACCESS_KEY_ID"
  plaintext_value = var.aws_access_key_id
}

resource "github_actions_secret" "aws-secret-access-key" {
  repository      = var.repo_name
  secret_name     = "AWS_SECRET_ACCESS_KEY"
  plaintext_value = var.aws_secret_access_key
}

resource "github_actions_secret" "gh-token" {
  repository      = var.repo_name
  secret_name     = "GH_TOKEN"
  plaintext_value = var.gh_token
}

resource "github_actions_secret" "repo-host" {
  repository      = var.repo_name
  secret_name     = "REPO_HOST"
  plaintext_value = var.repo_host
}

resource "github_actions_secret" "repo-user" {
  repository      = var.repo_name
  secret_name     = "REPO_USER"
  plaintext_value = var.repo_user
}

resource "github_actions_secret" "repo-key" {
  repository      = var.repo_name
  secret_name     = "REPO_KEY"
  plaintext_value = var.repo_key
}
