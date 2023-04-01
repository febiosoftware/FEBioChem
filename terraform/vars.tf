variable "group_name" {
  type    = string
  default = "github-ci-group"
}

variable "repo_name" {
  type    = string
  default = "FEBioChem"
}

variable "gh_token" {
  type      = string
  sensitive = true
}

variable "repo_host" {
  type      = string
  sensitive = false
}

variable "repo_user" {
  type      = string
  sensitive = false
}

variable "repo_key" {
  type      = string
  sensitive = true
}

variable "repo_action_permissions" {
  type    = list(string)
  default = ["aws/*", "tundranerd/*"]
}
