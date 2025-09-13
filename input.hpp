#pragma once
#include <memory>
#include <string_view>
#include <variant>
#include <vector>

enum class TokenType {
    Add,
    Sub,
    Mul,
    Pow,
    Var,
    Const,
    LPar,
    RPar,
};

struct Token {
    TokenType type;
    std::string_view value;
};

enum class ASTNodeType {
    Add,
    Sub,
    Mul,
    Pow,
    Var,
    Const,
};

struct ASTNode;
using ASTNodePtrType = std::unique_ptr<ASTNode>;
using ASTNodeForkType = std::pair<ASTNodePtrType, ASTNodePtrType>;
using ASTNodeConstType = double;
using ASTNodeVarType = char;

struct ASTNode {
    ASTNodeType type;
    std::variant<ASTNodeForkType, ASTNodeVarType, ASTNodeConstType> value;
};

/// Parse a vector of tokens into an abstract syntax tree.
std::unique_ptr<ASTNode> parse_tokens(std::vector<Token> tokens);

/// Tokenize a string expression into a vector of tokens.
std::vector<Token> tokenize(std::string_view expression);
