# Security Policy

## Supported Versions

| Version | Supported |
|---------|-----------|
| 2.0.x   | Yes       |
| < 2.0   | No        |

## Reporting a Vulnerability

If you discover a security vulnerability, please report it responsibly by
emailing <email@gdock.org>. Do **not** open a public issue.

You should expect an initial response within 7 days. The report will be
investigated and you will be kept informed of progress toward a fix.

## Scope

gdock is a command-line scientific tool that processes PDB files locally. It
does not handle authentication, network services, or sensitive user data.
Relevant security concerns include:

- Malformed PDB input causing crashes or unexpected behavior
- Path traversal via output directory options
- Dependency vulnerabilities
