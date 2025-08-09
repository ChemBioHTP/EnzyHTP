# Gemini PR Review Guide

This document outlines the process for Gemini to follow when reviewing pull requests for the EnzyHTP project. The goal is to ensure code quality, consistency, and adherence to project standards.

## 1. Understand the Context

Before diving into the code, it's crucial to understand the purpose and scope of the pull request.

1.  **Read the PR Description:** Start by thoroughly reading the pull request title, description, and any linked issues. This provides the context for the changes.
2.  **Identify Changed Files:** Get a list of all files modified in the PR. This will help you focus your review.
3.  **Review the Diff:** Use `git diff` or a similar tool to see the exact changes made to the code. Pay attention to not just what was added, but also what was removed or modified.

## 2. Code Review Checklist

Use the following checklist to guide your review of the code changes.

### Functionality and Correctness
- **Does it work?** Does the code achieve the stated goal of the PR?
- **Are there bugs?** Look for logical errors, edge cases that might not be handled, race conditions, or potential crashes.
- **Is it efficient?** Consider the performance implications of the changes. Is there a more efficient way to achieve the same result?

### Code Quality and Style
- **Readability:** Is the code clear, concise, and easy to understand? Are variable and function names descriptive?
- **Import Order:** Are the imports correctly organized according to the project conventions (standard library, third-party, then local application imports)?
- **Type Hinting:** Are functions properly type-hinted for clarity and static analysis?
- **Simplicity:** Is the code overly complex? Could it be simplified while still achieving the same goal?

### Architecture and Design
- **Adherence to Patterns:** Do the changes follow the established architectural patterns of the project (e.g., operating on `Structure` objects, using the interface layer for external tools)?
- **Modularity:** Is the code well-structured and modular? Does it have a clear separation of concerns?
- **Dependencies:** Does the PR introduce new external dependencies? If so, are they necessary and justified?

### Testing
- **Test Coverage:** Does the PR include new tests for the added or modified code?
- **Test Quality:** Are the tests well-written and meaningful? Do they cover both success and failure cases?
- **Passing Tests:** Do all existing and new tests pass? Run the relevant tests to verify. **Remember to run specific tests rather than the full suite.**

### Documentation
- **Docstrings:** Are new functions, classes, and modules documented with clear and informative docstrings?
- **Comments:** Are there comments where the code is complex or non-obvious? The comments should explain *why* the code is written a certain way, not *what* it does.

### Security
- **Input Validation:** Is all external input properly validated and sanitized?
- **Secrets:** Are any secrets or sensitive information being exposed?
- **Dependencies:** Are the dependencies secure and up-to-date?

## 3. Providing Feedback

When you provide feedback, follow these guidelines:

- **Be Constructive:** Your feedback should be helpful and aimed at improving the code.
- **Be Specific:** Reference the exact file and line number you are commenting on.
- **Explain Your Reasoning:** Clearly explain *why* you are suggesting a change.
- **Offer Suggestions:** When possible, provide concrete suggestions for improvement.
- **Prioritize:** Distinguish between critical issues that must be fixed and minor suggestions that are nice to have.

## 4. Review Workflow

1.  Start by stating the purpose of the PR.
2.  Use `git diff` to analyze the changes.
3.  Go through the checklist above for each changed file.
4.  Run the relevant tests to ensure that the changes are working as expected and do not break existing functionality.
5.  Summarize your findings in a clear and concise manner.
6.  If there are issues, provide specific feedback and suggestions for improvement.
