terraform {
  required_providers {
    github = {
      source  = "integrations/github"
      version = "~> 5.0"
    }
  }

  backend "s3" {
    bucket               = "febio-tf-state"
    key                  = "febio.tfstate"
    region               = "us-east-1"
    workspace_key_prefix = "febio"
    profile              = "febio"
    encrypt              = true
  }
}
