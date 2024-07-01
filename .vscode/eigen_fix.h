// file needed to prevent intellisense errors
// https://github.com/microsoft/vscode-cpptools/issues/7413
// https://github.com/microsoft/vscode-cpptools/issues/7413#:~:text=fixed%20after%20following-,%237413%20(comment),-%3F%20If%20not%2C%20can
#if __INTELLISENSE__
#undef __ARM_NEON
#undef __ARM_NEON__
#endif
