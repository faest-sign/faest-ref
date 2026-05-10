# Security Policy

## Reporting a Vulnerability

If you believe you have found a security issue, cryptographic correctness issue, side-channel concern, specification/implementation divergence, or vulnerability affecting this repository, please report it responsibly.

Please do not open a public issue with full exploit details or reproducer code before maintainers have had a chance to review the report.

Preferred disclosure contact:

```text
TODO: maintainers should add a private security contact email or enable GitHub private vulnerability reporting.
````

When reporting an issue, please include as much detail as possible:

* affected commit, branch, or release;
* affected parameter set, if applicable;
* affected files and functions;
* expected behavior according to the specification;
* observed implementation behavior;
* reproduction steps or proof-of-concept, if available;
* impact assessment;
* whether the issue appears to affect interoperability, correctness, side-channel resistance, or cryptographic security.

## Scope

Security reports may include, but are not limited to:

* signature verification bypasses;
* malformed signature acceptance;
* Fiat-Shamir transcript binding issues;
* parameter-set or domain-separation issues;
* specification-to-implementation divergences;
* non-canonical encoding or parsing ambiguity;
* randomness, nonce, salt, or seed-handling issues;
* side-channel concerns;
* memory safety issues;
* constant-time violations;
* test vector or KAT inconsistencies.

## Public Discussion

If a report does not require embargoed handling, maintainers may ask the reporter to open a public issue after initial review.

If the issue is confirmed and security-relevant, maintainers may coordinate a patch, advisory, credit, and disclosure timeline.

## Supported Versions

This repository currently tracks the reference implementation. Maintainers should update this section if specific releases or branches are supported.

| Version / Branch | Supported          |
| ---------------- | ------------------ |
| `main`           | Yes                |
| other branches   | Maintainer-defined |

## Disclosure Timeline

Maintainers should define the expected response and disclosure timeline.

Suggested default:

* acknowledge receipt within 7 days;
* provide an initial assessment within 30 days;
* coordinate public disclosure once a fix, specification clarification, or advisory is ready.
