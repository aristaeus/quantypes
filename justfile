test-default:
    cargo test

test-no-features:
    cargo test --no-default-features

test: test-default test-no-features
